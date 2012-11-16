function [] = ODE_Roberta_Simp_mMALA_PopMCMC(OutputID)

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
MCMC_Options.NumOfTemps            = 1;

% Set up temperature schedule
MCMC_Options.Temperatures          = 1; %(0:(1/(MCMC_Options.NumOfTemps-1)):1).^5; % Could use automatic temp schedule; see Friel et al. (2012).

% Set number of burn-in and posterior samples, and how often to save
MCMC_Options.NumOfBurnInSamples    = 100;
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
Model.Type                         = 'ODE';

% Name of the mathematical model to evaluate
Model.Name                         = 'SBRoberta';

% Choose sampling procedure for parameters
Model.Sampler                      = 'Simp_mMALA';

% Choose evaluation function to use - for calculating LL etc
Model.Func_Evaluate                = @ODE_Evaluate_Simp_mMALA_Numerical;
Model.Func_Proposal                = @ODE_Proposal_Simp_mMALA;

% Set up parameter names for reference
Model.ParaNames                    = {'', '', ''};
Model.NumOfParas                   = 80;

% Set up starting values for all temperatures
% Whether to use log space for parameter values - helps with very small/large values
Model.Log10Space                   = true;
Model.Paras                        = log10([0.250000000000000,0.0500000000000000,0.750000000000000,0.0500000000000000,0.750000000000000,0.0200000000000000,0.0100000000000000,0.0250000000000000,0.0250000000000000,0.0100000000000000,0.0250000000000000,10,0.0300000000000000,0.0300000000000000,0.00500000000000000,0.0500000000000000,0.0300000000000000,2.50000000000000,100,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.200000000000000,0.0200000000000000,0.550000000000000,0.0500000000000000,0.0100000000000000,0.750000000000000,0.750000000000000,0.0100000000000000,0.0250000000000000,0.750000000000000,10,0.0250000000000000,0.0100000000000000,0.0750000000000000,0.0500000000000000,0.0250000000000000,10,0.250000000000000,8,0.0250000000000000,10,9,0.250000000000000,8,0.0250000000000000,15,0.750000000000000,15,0.0250000000000000,15,0.750000000000000,15,0.0250000000000000,15,0.500000000000000,15,0.0250000000000000,15,0.500000000000000,15,0.0250000000000000,15,0.250000000000000,15,0.0250000000000000,15,0.0250000000000000,8,0.00200000000000000,0.0500000000000000,15,0.200000000000000,0.100000000000000,0.100000000000000,0.0500000000000000,0.0100000000000000,0.750000000000000;]);
%Model.Paras = [];
Model.StepSize                     = 1;

% Set up blocking for parameters - either fixed blocks or random blocks
% (This is used for higher-dimensional problems with ill-conditioned metric
% tensors)

% Either 'fixed':
%{
Model.Blocking.Type                = 'fixed';
Model.Blocking.NumOfBlocks         = 1;
Model.Blocking.Blocks{1}           = [1:80];
%}

% or 'random'
%
Model.Blocking.Type                = 'random'; % Either 'fixed' or 'random'
Model.Blocking.NumOfBlocks         = 20;
Model.Blocking.Blocks              = [ones(1, 20)*4]; % Must sum to number of parameters
%}

% Set up priors for each of the parameters
if Model.Log10Space
    
    % Set up priors for log space
    
    for i = 1:Model.NumOfParas
        Model.Func_Evaluate_Prior{i} = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{i}   = @Sample_Prior_Uniform;
        Model.Priors{i}.Paras        = [-4 3];
    end
    
else
    
    % Set up priors for standard space
    
    for i = 1:Model.NumOfParas
        Model.Func_Evaluate_Prior{i} = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{i}   = @Sample_Prior_Uniform;
        Model.Priors{i}.Paras        = [0 200];
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE model settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up names of ODE solvers for 1st order sensitivities
Model.SBModelName_1st_Paras        = '';

% Set up starting values for all temperatures
% Use log space for initial conditions - helps with very small/large values
Model.Log10Space_ICs               = false;
Model.ICs                          = [1000,0,500,0,50,0,300,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,50,0,50,0,50,0,50,0,0,100,0,0,100,0,0,0,0,0,50,0;];

% Set total number of species in ODE system
Model.NumOfSpecies                 = 42;

% Set which species are observed and unobserved
Model.ObservedSpecies              = 1:42;
Model.UnobservedSpecies            = [];

% Set solver tolerances
Model.Options.abstol               = 1e-8;
Model.Options.reltol               = 1e-8;

% Decide whether to infer initial conditions
Model.InferICs                     = false; % otherwise keep them fixed at above values

% Set up blocking for ODE initial conditions - either fixed blocks or random blocks
% (This is used for higher-dimensional problems with ill-conditioned metric
% tensors)
if Model.InferICs
    
    % Choose sampling procedure for initial conditions
    Model.Sampler_ICs                  = 'MH';
    
    % Choose evaluation function to use - for calculating LL etc
    Model.Func_Evaluate_ICs            = @ODE_Evaluate_ICs_MH;
    Model.Func_Proposal_ICs            = @ODE_Proposal_ICs_MH;
    
    Model.Upper_ICs                    = 0.5;
    Model.Lower_ICs                    = 0.2;
    
    % Set up names of ODE solvers for 1st order sensitivities for ICs
    Model.SBModelName_1st_ICs          = '';
    
    Model.StepSize_ICs                 = 0.05;
    
    % Either 'fixed':
    %
    Model.Blocking_ICs.Type            = 'fixed';
    Model.Blocking_ICs.NumOfBlocks     = 1;
    Model.Blocking_ICs.Blocks{1}       = [1 2];
    %}
    
    % or 'random'
    %{
    Model.Blocking_ICs.Type            = 'random'; % Either 'fixed' or 'random'
    Model.Blocking_ICs.NumOfBlocks     = 2;
    Model.Blocking_ICs.Blocks          = [1 1]; % Must sum to number of initial conditions
    %}
    
    % Set up priors for initial conditions
    if Model.Log10Space_ICs
        
        % Set up priors for log space
        
        for i = 1:Model.NumOfSpecies
            Model.Func_Evaluate_Prior_ICs{i}   = @Evaluate_Prior_Normal;
            Model.Func_Sample_Prior_ICs{i}     = @Sample_Prior_Normal;
            Model.Priors_ICs{i}.Paras          = [0 1];
        end
        
    else
        
        % Set up priors for standard space
        
        for i = 1:Model.NumOfSpecies
            Model.Func_Evaluate_Prior_ICs{i}   = @Evaluate_Prior_Uniform;
            Model.Func_Sample_Prior_ICs{i}     = @Sample_Prior_Uniform;
            Model.Priors_ICs{i}.Paras          = [-5 5];
        end
        
    end
    
    
end


% Set up the data for the ODE model

% Open the given file
InputID                            = 'ODE_Roberta_FullyObs_5pc.mat';
InputData                          = open(InputID);

% Save the data and timepoints
Model.Data                         = InputData.Data;
Model.TimePoints                   = InputData.TimePoints;
NoiseVariance                      = InputData.NoiseVariance;
Model.NumOfTimePoints              = size(Model.TimePoints, 1);

% Save the covariance matrices for calculating the likelihood
for i = Model.ObservedSpecies
    Model.NoiseCov{i} = eye(Model.NumOfTimePoints)*NoiseVariance(i);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call main population MCMC routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Population_MCMC(MCMC_Options, Model);



end