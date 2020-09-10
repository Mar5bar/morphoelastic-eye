function H = scleralReferenceThickness(Sigma,referenceRadius,option)
% scleralReferenceThickness	Compute the reference thickness of the sclera.
% 
% 	H = scleralReferenceThickness(Sigma,R,Option) computes the reference
% 	arclength at Sigma for a spherical reference configuration of radius R.
% 
% 	Option is one of 'uniform' and 'non-uniform', generating the appropriate
% 	reference uniform or non-uniform thickness.
	
	% For brevity.
	R = referenceRadius;

	% Return the reference thickness at the Sigmas.
	% This may be uniform or non-uniform.
	switch option
	case 'uniform'
		% -------
		% Uniform
		% -------
		H = 0*Sigma + 0.65;

	case 'non-uniform'
		% -----------
		% Non-uniform
		% -----------
		H = (Sigma < R*pi/2) .* (0.65  - tanh((Sigma - pi*R/6 )/1.5) * 0.2) + ...
	    	(Sigma > R*pi/2) .* (0.475 + tanh(Sigma - 0.7*pi*R) * 0.025);
	otherwise
		error('Invalid reference thickness option. Option must be uniform or non-uniform.')
    end
end