function [ Chain, Success ] = ODE_Evaluate_Simp_mMALA_Obs_Numerical( Model, Chain )

% Inputs: Model object, Chain object
% Output: Chain object, Success flag


% Initialise success flag
Success = 1;


NumOfParas   = Model.NumOfParas; % Number of parameters
NumOfSpecies = Model.ModelSpecific.NumOfSpecies; % Number of species



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


if Model.ModelSpecific.Log10Space_ICs
    ICs = 10.^ICs;
end

% Assume all time points are the same, for now...
TimePoints = Model.TimePoints(:, 1);

try
    % Solve the ODEs
    SimData = feval(Model.Name, TimePoints, ICs, Paras, Model.ModelSpecific.Options);
catch
    Success = 0;
    return
end



%%% Calculate the log-likelihood %%%

% Get the simulated states
StateEstimates = SimData.reactionvalues;

NumOfReactions = size(StateEstimates,2);

% Calculate Log-Likelihood for each species
for n = 1:NumOfReactions
    LL(n) = Normal_LogPDF(StateEstimates(:,n), Model.Data(:,n), Model.NoiseCov{n});
end

Chain.LL = sum(LL);



%%% Calculate the gradient of the log-likelihood %%%

Epsilon = 1e-6;

% Check tol is low enough
if (Epsilon - Model.ModelSpecific.Options.abstol) <= 0
    disp('Error ODE solver needs smaller tolerance for stable numerical estimates.')
    Success = 0;
    return
end

% Numerically estimate the required ODE sensitivities
for ParaNum = 1:Chain.CurrentBlockSize
    
    CurrentPara                = Chain.CurrentBlock(ParaNum);
    
    NewParas                   = Chain.Paras;
    NewParas(CurrentPara)      = NewParas(CurrentPara) + Epsilon;
    
    
    if Model.ModelSpecific.InferICs
        % Separate parameters and initial conditions
        NewICs   = NewParas((Model.NumOfParas-Model.ModelSpecific.NumOfSpecies+1):end);
        NewParas = NewParas(1:(Model.NumOfParas-Model.ModelSpecific.NumOfSpecies));
        
    else
        % Just parameters so set up ICs
        NewICs   = Model.ModelSpecific.ICs;
    end
    
    
    NewSimData                 = feval(Model.Name, TimePoints, NewICs, NewParas, Model.ModelSpecific.Options);
    
    Sensitivities{CurrentPara} = (NewSimData.reactionvalues - StateEstimates)./Epsilon;
    
end



% Calculate gradients for each of the parameters i.e. d(LL)/d(Parameter)
Chain.GradLL       = zeros(1, Chain.CurrentBlockSize);
Chain.GradLogPrior = zeros(1, Chain.CurrentBlockSize);


for ParaNum = 1:Chain.CurrentBlockSize
    
    CurrentPara = Chain.CurrentBlock(ParaNum);
    
    for i = Model.ModelSpecific.ObservedSpecies
        Chain.GradLL(ParaNum) = Chain.GradLL(ParaNum) - ((StateEstimates(:,i)-Model.Data(:,i))'/Model.NoiseCov{i})*Sensitivities{CurrentPara}(:,i);
    end
    
    % Add log transformation if necessary
    if Model.Log10Space
        Chain.GradLL(ParaNum) = Chain.GradLL(ParaNum)*( 10.^Chain.Paras(CurrentPara)*log(10) );
    end
    
    % Calculate the gradient of the log prior
    Chain.GradLogPrior(ParaNum) = Model.Func_Evaluate_Prior{CurrentPara}(Model.Priors{CurrentPara}.Paras, Chain.Paras(CurrentPara), 1);
end



%%% Calculate metric tensor %%%

Chain.FI = zeros(Chain.CurrentBlockSize);

for i = 1:Chain.CurrentBlockSize
    for j = i:Chain.CurrentBlockSize
        
        CurrentPara_i = Chain.CurrentBlock(i);
        CurrentPara_j = Chain.CurrentBlock(j);
        
        for SpeciesNum = Model.ModelSpecific.ObservedSpecies
            Chain.FI(i,j) = Chain.FI(i,j) + ((Sensitivities{CurrentPara_i}(:,SpeciesNum)'/Model.NoiseCov{SpeciesNum})*Sensitivities{CurrentPara_j}(:,SpeciesNum));
        end
        
        % Add log transformation if necessary
        if Model.Log10Space
            Chain.FI(i,j) = Chain.FI(i,j)*( 10.^Chain.Paras(CurrentPara_i)*log(10) )*( 10.^Chain.Paras(CurrentPara_j)*log(10) );
        end
        
        Chain.FI(j,i) = Chain.FI(i,j);
    end
end


Chain.HessianLogPrior = zeros(1, Chain.CurrentBlockSize);

% Calculate the Hessian of the log prior
for ParaNum = 1:Chain.CurrentBlockSize
    
    CurrentPara                    = Chain.CurrentBlock(ParaNum);
    
    Chain.HessianLogPrior(ParaNum) = Model.Func_Evaluate_Prior{CurrentPara}(Model.Priors{CurrentPara}.Paras, Chain.Paras(CurrentPara), 2);
end





end

