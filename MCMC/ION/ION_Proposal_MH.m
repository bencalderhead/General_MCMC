function [ Proposal_Mean, Proposal_Covariance ] = ION_Proposal_MH( Model, Chain )


Proposal_Mean    = Chain.Paras(Chain.CurrentBlock);

% Either
Proposal_Covariance = eye(Chain.CurrentBlockSize)*(Chain.StepSize(Chain.CurrentBlockNum)^2);
% or
%Proposal_Covariance = CholG\(CholG'\( diag(ones(1, Chain.NumOfParas)*(Chain.StepSize^2)) )); % Stepsize is squared in the covariance term

end

