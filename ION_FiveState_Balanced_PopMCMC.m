function [] = ION_FiveState_Balanced_PopMCMC(OutputID)

% Add all sub-folders for Matlab
addpath(genpath('./'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up MCMC options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set output file for saving samples
MCMC_Options.OutputID              = OutputID;

% Decide whether to save burn-in
MCMC_Options.SaveBurnIn            = true;

% Set number of tempered distributions to use
MCMC_Options.NumOfTemps            = 20;

% Set up temperature schedule
MCMC_Options.Temperatures          = (0:(1/(MCMC_Options.NumOfTemps-1)):1).^5; % Could use automatic temp schedule; see Friel et al. (2012).

% Set number of burn-in and posterior samples, and how often to save
MCMC_Options.NumOfBurnInSamples    = 2000;
MCMC_Options.NumOfPosteriorSamples = 5000;
MCMC_Options.PosteriorSaveSize     = 5000; % This is changed from .SaveFullEvery

% Set iteration interval for adapting stepsizes
MCMC_Options.AdaptRate             = 50;
MCMC_Options.Upper                 = 0.5;
MCMC_Options.Lower                 = 0.2;

% Set number of cores to use
MCMC_Options.NumOfCores            = 2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common model settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose model type - either ODE or ION
Model.Type                         = 'ION';

% Name of the mathematical model to evaluate
Model.Name                         = 'Five_State_Balanced';

% Choose sampling procedure for parameters
Model.Sampler                      = 'MH';

% Choose evaluation function to use - for calculating LL etc
Model.Func_Evaluate                = @ION_Evaluate_MH;
Model.Func_Proposal                = @ION_Proposal_MH;

% Set up parameter names for reference
Model.ParaNames                    = {'K_+1', 'K_-1', 'K_+2', 'K_-2', 'Beta_1', 'Alpha_1', 'Beta_2', 'Alpha_2', 'K*_+2'};
Model.NumOfParas                   = 9;

% Set up starting values for all temperatures
% Whether to use log space for parameter values - helps with very small/large values
Model.Log10Space                   = true;
Model.Paras                        = log10([5e7, 2000, 5e8, 2000, 15, 3000, 15000, 500, 5e8]);
Model.StepSize                     = 0.01;

% Set up blocking for parameters - either fixed blocks or random blocks
% (This is used for higher-dimensional problems with ill-conditioned metric
% tensors)

% Either 'fixed':
%
Model.Blocking.Type                = 'fixed';
Model.Blocking.NumOfBlocks         = 1;
Model.Blocking.Blocks{1}           = [1:9];
%}

% or 'random'
%{
Model.Blocking.Type                = 'random'; % Either 'fixed' or 'random'
Model.Blocking.NumOfBlocks         = 3;
Model.Blocking.Blocks              = [1 1 1]; % Must sum to number of parameters
%}

% Set up priors for each of the parameters
if Model.Log10Space
    
    % Set up priors for log space
    
    for i = 1:length(Model.Paras)
        Model.Func_Evaluate_Prior{i} = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{i}   = @Sample_Prior_Uniform;
        Model.Priors{i}.Paras        = [-2 10];
    end
    
else
    
    % Set up priors for standard space
    
    for i = 1:length(Model.Paras)
        Model.Func_Evaluate_Prior{i} = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{i}   = @Sample_Prior_Uniform;
        Model.Priors{i}.Paras        = [1e-2 1e10];
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ION model settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set total number of species in ODE system
Model.ModelSpecific.NumOfStates                  = 5;

% Set which species are observed and unobserved
Model.ModelSpecific.ClosedStates                 = [3 4 5];
Model.ModelSpecific.OpenStates                   = [1 2];


%%% Initialise any chain specific values %%%
Model.ChainSpecific = [];


% Set up the data for the ODE model

% Open the given file
InputID                            = 'ION_Five_State_Balanced_Data_Long.mat';
InputData                          = open(InputID);

% Save the data and timepoints
Model.Data                         = InputData.Data;
Model.TimePoints                   = InputData.TimePoints;
Model.NumOfTimePoints              = size(Model.TimePoints, 1);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call main population MCMC routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Population_MCMC(MCMC_Options, Model);



end