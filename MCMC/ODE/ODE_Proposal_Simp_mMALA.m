function [ Proposal_Mean, Proposal_Covariance ] = ODE_Proposal_Simp_mMALA( Model, Chain )

StepSize     = Chain.StepSize(Chain.CurrentBlockNum);
CurrentParas = Chain.CurrentBlock;
NumOfParas   = Chain.CurrentBlockSize;

% Calculate posterior gradient and metric tensor
Posterior_Grad   = Chain.GradLL*Chain.Temp + Chain.GradLogPrior;
Posterior_G      = Chain.FI*Chain.Temp - diag(Chain.HessianLogPrior);


% Calculate cholesky
CholG = chol(Posterior_G);
                
NaturalGradient  = (CholG\(CholG'\Posterior_Grad'))';
Proposal_Mean    = Chain.Paras(CurrentParas) + (StepSize^2/2)*NaturalGradient;

% Either
%Proposal_Covariance = inv(Posterior_G)*(Chain.StepSize^2);
% or
Proposal_Covariance = CholG\(CholG'\( diag(ones(1, NumOfParas)*(StepSize^2)) )); % Stepsize is squared in the covariance term

end

