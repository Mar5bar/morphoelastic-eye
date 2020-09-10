function gc = intrinsicGrowthCapacity(Sigma,Gamma,delta,eta0)
% intrinsicGrowthCapacity Compute the intrinsic growth capacity of the shell.
% 
% 	gc = intrinsicGrowthCapacity(Sigma,Gamma,delta,eta0) computes the
% 	intrinsic growth capacity of the shell at reference arclengths Sigma for a
% 	region at location Gamma of extent delta with scaling eta0.

	gc = eta0 * (tanh((Gamma - Sigma) / delta) + 1) / 2;

end