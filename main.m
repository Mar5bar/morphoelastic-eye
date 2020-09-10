% This code accompanies the manuscript 
%   `A shell model of eye growth and elasticity'
% by 

% This script will simulate the growth of a morphoelastic eye-like shell in
% response to optical feedback, as described in detail in the accompanying
% publication. Parameter naming convention is as in the accompanying
% manuscript, or more verbose to aid readability or if ambiguous. This and all
% related code is available under a CC BY 4.0 Attribution License.

% -----------------
% Solver parameters
% -----------------
% Discretisation level in arclength.
n_arclength=501;

% Discretisation level in time.
n_time=200;

% Level of truncation in arclength around s = sigma = Sigma = 0 to avoid
% numerical singularity, with the domain truncated at sigma = epsilon.
epsilon = 1e-3;

% Solver tolerance in the quasi-steady boundary value problem.
bvp_tolerance = 1e-6;

% -------------------------------------
% Shape and material parameters
% -------------------------------------
% Proportion of hemispherical initial reference configuration occupied by the
% sclera, at which point it smoothly attaches to the cornea. Between 0 and 1.
scleral_reference_proportion = 0.8;

% Radius of the sclera in the initial reference configuration. Note that this
% is shared by the reference configuration of the cornea.
scleral_reference_radius = 10;

% From these values, compute the reference arclength, i.e. the upper bound for
% Sigma.
L = pi * scleral_reference_proportion * scleral_reference_radius;

% Interocular pressure, constant throughout.
IOP = 2;

% Tensile strength of the sclera.
C = 100;

% Tensile strength of the cornea.
cornea_C = 100;

% Uniform thickness of the cornea in its reference configuration.
corneal_reference_thickness = 0.65;

% Strength of fibre reinforcement in the sclera.
D = 0;

% Bending stiffness of the sclera.
EB = 18;

% Uniform or non-uniform reference scleral thickness? Details can be altered
% in scleralReferenceThickness.m
% thicknessOption = 'uniform';
thicknessOption = 'non-uniform';

% -----------------
% Growth parameters
% -----------------
% Location of the region of the sclera that can grow.
Gamma = 8;

% Spatial extent of the region of the sclera that can grow.
delta = 4;

% Rate of scleral growth scaling.
eta0 = 70;

% -------------
% Initial setup
% -------------

% Compute the timestep for later use.
timestep = 1/n_time;

% Z coordinate at which the undeformed sclera attaches to the undeformed
% cornea, required in order to correctly set the origin in the grown but
% undeformed configuration.
corneal_reference_z = scleral_reference_radius*(1+cos(scleral_reference_proportion*pi));

% Compute the constant corneal stretch.
const = IOP * scleral_reference_radius / (4 * corneal_reference_thickness * cornea_C);
% Solve to high accuracy between 1 and 7^(1/6) for the stretch.
options = optimset('TolFun',1e-10,'TolX',1e-10);
corneal_stretch_bounds = [1, 7^(1/6)];
corneal_stretch = fzero(@(lambda) 1./lambda - 1./lambda.^7 - const, corneal_stretch_bounds);

% Compute the coordinates of the endpoint of the deformed cornea, to be used
% as boundary conditions for the sclera.
corneal_theta = scleral_reference_proportion * pi;
corneal_r = corneal_stretch * scleral_reference_radius * sin(scleral_reference_proportion * pi);
corneal_z = corneal_stretch * scleral_reference_radius * (1 + cos(scleral_reference_proportion * pi));

% Discretise the reference arclength Sigma, which runs between epsilon and L.
Sigma = linspace(epsilon, L, n_arclength);

% Compute the R coordinates for the scleral reference configuration.
R = scleral_reference_radius * sin(Sigma / scleral_reference_radius);

% Get the reference scleral thickness.
H = scleralReferenceThickness(Sigma, scleral_reference_radius, thicknessOption);

% Set up fibre orientation, psi.
psi = fibreOrientation(Sigma,R); 

% Set up the intrinsic growth capacity.
gc = intrinsicGrowthCapacity(Sigma, Gamma, delta, eta0);

% ----------------------------------------------------
% Compute the initial ungrown, deformed configuration.
% ----------------------------------------------------

% Initially, before growth has occured, sigma = Sigma and rho = R.
sigma = Sigma;
rho = R;

% The solver requires an initial guess for all quantities.
init_guess = struct();
% The bvp will have sigma as its independent variable.
init_guess.x = sigma;
% There will be five unknown dependent variables to solve for.
init_guess.y = zeros(5,length(sigma));
% These are ordered as kappa, alpha_s, r, theta, and Q.
init_guess.y(1,:) = ones(length(sigma),1) / scleral_reference_radius;
init_guess.y(2,:) = ones(length(sigma),1);
init_guess.y(3,:) = scleral_reference_radius * sin(sigma / scleral_reference_radius);
init_guess.y(4,:) = sigma / scleral_reference_radius;
init_guess.y(5,:) = 0;

% Solve the elastic bvp. The solution will be evaluated at the material points
% specified by the discretisation of Sigma, stored in sigma.
sol = StretchBvp(sigma,rho,init_guess,corneal_r,corneal_theta,IOP,C,H,psi,D,EB,bvp_tolerance);

% Extract the solution components kappa, alpha_s, r, theta, and Q. Note that
% these are evaluated at the material points sigma(Sigma).
kappa = sol(1,:);
alpha_s = sol(2,:);
r = sol(3,:);
theta = sol(4,:);
Q = sol(5,:);

% Compute the elastic stretch in the phi and normal direction.
alpha_phi = r ./ rho;
alpha_n = 1 ./ (alpha_s .* alpha_phi);

% Recover the z coordinate from the deformed solution, setting the origin at
% the anterior eye with the z coordinate of the corneal attachment.
z = cumtrapz(sigma, - alpha_s .* sin(theta));
z = z - z(end) + corneal_z;

% Recover s, the grown deformed arclength, at these material points.
s = cumtrapz(sigma,alpha_s) + alpha_s(1) * epsilon;

% Compute the deformed scleral thickness, h = H * alpha_n.
h = H .* alpha_n;

% ----------------------------------------------------------
% Compute the growth capacity of this deformed configuration
% ----------------------------------------------------------

% Load in the information needed to construct the target surface, the
% best-focus surface for blue light. This is defined by a number of
% quantities: BFS_P, BFS_point, and BFS_max_angle, where BFS_point is a
% reference point (0,0,z) about which the surface has a polar form with
% parameters by BFS_P, and BFS_max_angle is the maximum angle relative to
% BFS_point up to which the surface is defined.
load('target_surface.mat')

% Compute the location of the retina, at the front face of the deformed
% sclera.
z_retina = z - h .* cos(theta) / 2;
r_retina = r - h .* sin(theta) / 2;

% Parameterise the surface of the retina by an angle about the BFS_point.
retina_angle = atan2(r_retina, z_retina - BFS_point);

% We can only measure the distance between the target surface and the retina
% where both are defined, with the limiting factor being the target surface.
mask = retina_angle < BFS_max_angle;
retina_angle = retina_angle(mask);
r_retina = r_retina(mask);
z_retina = z_retina(mask);

% Compute the radius of the retina relative to the point (0,0,BFS_point).
radius_retina = sqrt(r_retina.^2 + (z_retina - BFS_point).^2);

% Reconstruct the target surface at the retina_angles using its quartic
% representation.
radius_BFS = BFS_P(1) * retina_angle.^4 + BFS_P(2) * retina_angle.^2 + BFS_P(3);

% Compute the distance between the two surfaces at each retina_angle-defined
% material point. This is a simple subtraction of their radii.
distance = radius_BFS - radius_retina;

% If the retina is behind the target surface, do not grow.
distance(distance < 0) = 0;

% Extrapolate beyond the target surface, taking the last value.
gv = zeros(1,length(sigma));
gv(mask) = distance;
gv(~mask) = distance(end);

% Form the growth rate, eta.
eta = gc .* gv;

% Save the configurations to file.
time_index = 0;
time = time_index * timestep;
save(['output_timestep_',num2str(time_index),'.mat'],'h','s','z','r','sigma','rho','eta','time')

% -----------------------
% Growth-deformation loop
% -----------------------


for time_index = 1 : n_time

    time = time_index * timestep;

    % Grow the configuration.
    rho = rho + timestep * eta .* max(r-rho,0);
    sigma = sigma + timestep * (cumtrapz(sigma, max(alpha_s - 1, 0) .* eta) + sigma(1) * eta(1) * max(alpha_s(1) - 1,0));
    
    % Update the initial guess for the upcoming bvp solver to be the previous solution, but with updated sigma.
    init_guess.x = sigma;
    init_guess.y = [kappa; alpha_s; r; theta; Q];

    % Solve the elastic bvp. The solution will again be evaluated at the
    % material points specified by the discretisation of Sigma, stored in
    % sigma.
    sol = StretchBvp(sigma,rho,init_guess,corneal_r,corneal_theta,IOP,C,H,psi,D,EB,bvp_tolerance);

    % Extract the solutions.
    kappa = sol(1,:);
    alpha_s = sol(2,:);
    r = sol(3,:);
    theta = sol(4,:);
    Q = sol(5,:);

    % Compute the elastic stretch in the phi and normal direction.
    alpha_phi = r ./ rho;
    alpha_n = 1 ./ (alpha_s .* alpha_phi);

    % Recover the z coordinate from the deformed solution, setting the origin at
    % the anterior eye with the z coordinate of the corneal attachment.
    z = cumtrapz(sigma, - alpha_s .* sin(theta));
    z = z - z(end) + corneal_z;

    % Recover s, the grown deformed arclength, at these material points.
    s = cumtrapz(sigma,alpha_s) + alpha_s(1) * epsilon;

    % Compute the deformed scleral thickness, h = H * alpha_n.
    h = H .* alpha_n;
    
    % --------------------------
    % Update the growth function
    % --------------------------

    % Compute the location of the retina, at the front face of the deformed
    % sclera.
    z_retina = z - h .* cos(theta) / 2;
    r_retina = r - h .* sin(theta) / 2;

    % Parameterise the surface of the retina by an angle about the BFS_point.
    retina_angle = atan2(r_retina, z_retina - BFS_point);

    % We can only measure the distance between the target surface and the retina
    % where both are defined, with the limiting factor being the target surface.
    mask = retina_angle < BFS_max_angle;
    retina_angle = retina_angle(mask);
    r_retina = r_retina(mask);
    z_retina = z_retina(mask);

    % Compute the radius of the retina relative to the point (0,0,BFS_point).
    radius_retina = sqrt(r_retina.^2 + (z_retina - BFS_point).^2);

    % Reconstruct the target surface at the retina_angles using its quartic
    % representation.
    radius_BFS = BFS_P(1) * retina_angle.^4 + BFS_P(2) * retina_angle.^2 + BFS_P(3);

    % Compute the distance between the two surfaces at each
    % retina_angle-defined material point. This is a simple subtraction of
    % their radii.
    distance = radius_BFS - radius_retina;

    % If the retina is behind the target surface, do not grow.
    distance(distance < 0) = 0;

    % Extrapolate beyond the target surface, taking the last value.
    gv = zeros(1,length(sigma));
    gv(mask) = distance;
    gv(~mask) = distance(end);

    % Form the growth rate, eta.
    eta = gc .* gv;

    % Save the configurations to file.
    save(['output_timestep_',num2str(time_index),'.mat'],'h','s','z','r','sigma','rho','eta','time')
 
end