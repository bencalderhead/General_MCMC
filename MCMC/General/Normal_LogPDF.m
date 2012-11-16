function [ LogResult ] = Normal_LogPDF( Mean, Mu, Covar )

% Calculate probability of Mean using *normalised* multivariate Gaussian

% Mean and Mu are column vectors (k by 1)
% Covar is (pos-def) covariance matrix (k by k)

k = length(Mean);

Diff = (Mean-Mu);

LogResult = -(k/2)*log(2*pi) - sum(log(diag(chol(Covar)))) -0.5*(( Diff'/Covar )*Diff);


end

