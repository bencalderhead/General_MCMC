function [ Chain, Success ] = ODE_Evaluate_MH( Model, Chain )

% Inputs: Model object, Chain object
% Output: Chain object, Success flag


% Initialise success flag
Success    = 1;


NumOfParas = Model.NumOfParas; % Number of parameters



%%% Calculate the prior for each parameter %%%

for ParaNum = 1:NumOfParas
    Prior(ParaNum) = Model.Func_Evaluate_Prior{ParaNum}(Model.Priors{ParaNum}.Paras, Chain.Paras(ParaNum), 0);
end

if min(Prior) > 0
    
    % Convert into log value
    Chain.LogPrior = log(Prior);
else
    
    % Proposed parameter values lie outside prior, so fail and exit
    Success = 0;
    return
end



%%% Solve augmented ODE to get 1st order sensitivities %%%

if Model.ModelSpecific.InferICs
    % Separate parameters and initial conditions
    Paras = Chain.Paras(1:(Model.NumOfParas-Model.ModelSpecific.NumOfSpecies));
    ICs   = Chain.Paras((Model.NumOfParas-Model.ModelSpecific.NumOfSpecies+1):end);
else
    % Just parameters
    Paras = Chain.Paras;
    ICs   = Model.ModelSpecific.ICs;
end
    

% Reparameterise if necessary
if Model.Log10Space
    Paras = 10.^Paras;
end

if Model.ModelSpecific.InferICs
    if Model.ModelSpecific.Log10Space_ICs
        ICs = 10.^ICs;
    end
end

% Assume all time points are the same, for now...
TimePoints = Model.TimePoints(:, 1);

try
    % Solve the 1st order equations
    SimData = feval(Model.Name, TimePoints, ICs, Paras, Model.ModelSpecific.Options);
catch
    Success = 0;
    Chain.LL = -1e300;
    return
end



%%% Calculate the log-likelihood %%%

% Get the simulated states
StateEstimates = SimData.statevalues;

% Calculate Log-Likelihood for each species
for n = Model.ModelSpecific.ObservedSpecies
    LL(n) = Normal_LogPDF(StateEstimates(:,n), Model.Data(:,n), Model.NoiseCov{n});
end

Chain.LL = sum(LL);





end

