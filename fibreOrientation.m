function psi = fibreOrientation(Sigma,scleral_R)
% fibreOrientation Compute the constant fibre orientation.
% 
% 	psi = fibreOrientation(Sigma,R) computes the fibre orientation at reference
% 	arclengths Sigma with accompanying radius R.

	% Find the equator of the reference configuration, where scleral_R is
	% equal to its maximum.
	[~,indexEquator] = max(scleral_R);
	SigmaEquator = Sigma(indexEquator);
	% Define intermediate quantities.
	b = SigmaEquator / 4;
	a = b / 5;
	% Compute psi as in the accompanying manuscript.
	psi = 0.25 + 0.22 * (Sigma >= SigmaEquator) .* cos(pi * (Sigma - SigmaEquator) / (Sigma(end) - SigmaEquator)) + ...
	     (Sigma >= b) .* (Sigma < SigmaEquator) .* (0.05 - 0.17 * cos(pi * (Sigma - b) / (SigmaEquator - b))) + ...
	     (Sigma >= a) .* (Sigma < b) .* (0.06 * cos(pi * (Sigma - a) / (b - a)) - 0.06);
	psi = pi * psi;
end