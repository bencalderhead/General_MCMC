function [ Proposal_Mean, Proposal_Covariance ] = ODE_Proposal_MH( Model, Chain )


Proposal_Mean       = Chain.Paras(Chain.CurrentBlock);

% Either
Proposal_Covariance = eye(Chain.CurrentBlockSize)*(Chain.StepSize(Chain.CurrentBlockNum)^2);

end

