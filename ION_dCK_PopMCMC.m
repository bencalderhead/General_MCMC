function [] = ION_dCK_PopMCMC(OutputID)

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
MCMC_Options.NumOfTemps            = 10;

% Set up temperature schedule
MCMC_Options.Temperatures          = (0:(1/(MCMC_Options.NumOfTemps-1)):1).^5; % Could use automatic temp schedule; see Friel et al. (2012).

% Set number of burn-in and posterior samples, and how often to save
MCMC_Options.NumOfBurnInSamples    = 1000;
MCMC_Options.NumOfPosteriorSamples = 1000;
MCMC_Options.PosteriorSaveSize     = 1000; % This is changed from .SaveFullEvery

% Set iteration interval for adapting stepsizes
MCMC_Options.AdaptRate             = 50;
MCMC_Options.Upper                 = 0.5;
MCMC_Options.Lower                 = 0.2;

% Set number of cores to use
MCMC_Options.NumOfCores            = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common model settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose model type - either ODE or ION
Model.Type                         = 'ION';

% Name of the mathematical model to evaluate
Model.Name                         = 'Castillo_Katz';

% Choose sampling procedure for parameters
Model.Sampler                      = 'MH';

% Choose evaluation function to use - for calculating LL etc
Model.Func_Evaluate                = @ION_Evaluate_MH;
Model.Func_Proposal                = @ION_Proposal_MH;

% Set up parameter names for reference
Model.ParaNames                    = {'K_1', 'K_2', 'Beta', 'Alpha'};
Model.NumOfParas                   = 4;

% Set up starting values for all temperatures
% Whether to use log space for parameter values - helps with very small/large values
Model.Log10Space                   = true;
Model.Paras                        = log10([100 100 1000 1000]);
Model.StepSize                     = 1;

% Set up blocking for parameters - either fixed blocks or random blocks
% (This is used for higher-dimensional problems with ill-conditioned metric
% tensors)

% Either 'fixed':
%
Model.Blocking.Type                = 'fixed';
Model.Blocking.NumOfBlocks         = 1;
Model.Blocking.Blocks{1}           = [1 2 3 4];
%Model.Blocking.Blocks{2}           = [2];
%Model.Blocking.Blocks{3}           = [3];
%Model.Blocking.Blocks{4}           = [4];
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
        Model.Priors{i}.Paras        = [-2 4];
    end
    
else
    
    % Set up priors for standard space
    
    for i = 1:length(Model.Paras)
        Model.Func_Evaluate_Prior{i} = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{i}   = @Sample_Prior_Uniform;
        Model.Priors{i}.Paras        = [0 10000];
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ION model settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set total number of species in ODE system
Model.ModelSpecific.NumOfStates                  = 3;

% Set which species are observed and unobserved
Model.ModelSpecific.ClosedStates                 = [1 2];
Model.ModelSpecific.OpenStates                   = [3];


%%% Initialise any chain specific values %%%
Model.ChainSpecific = [];


% Set up the data for the ION model

% Open the given file
InputID                            = 'ION_dCK_0,5s.mat';
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
