function [] = Locke_Benchmark_1_Simp_mMALA(OutputID)

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
MCMC_Options.NumOfTemps            = 1; %20;

% Set up temperature schedule
MCMC_Options.Temperatures          = 1; %(0:(1/(MCMC_Options.NumOfTemps-1)):1).^5; % Could use automatic temp schedule; see Friel et al. (2012).

% Set number of burn-in and posterior samples, and how often to save
MCMC_Options.NumOfBurnInSamples    = 5000;
MCMC_Options.NumOfPosteriorSamples = 20000;
MCMC_Options.PosteriorSaveSize     = 20000; % This is changed from .SaveFullEvery

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
Model.Name                         = 'SBLocke_2005a_Full';

% Choose sampling procedure for parameters
Model.Sampler                      = 'Simp_mMALA';

% Choose evaluation function to use - for calculating LL etc
Model.Func_Evaluate                = @ODE_Evaluate_Simp_mMALA_Numerical;
Model.Func_Proposal                = @ODE_Proposal_Simp_mMALA;

% Set up parameter names for reference
Model.ParaNames                    = {'g1', 'g2', 'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'n1', 'n2', 'p1', 'p2', 'r1', 'r2', 'r3', 'r4'};
Model.NumOfParas                   = 22;

% Set up starting values for all temperatures
% Whether to use log space for parameter values - helps with very small/large values
Model.Log10Space                   = false;
Model.Paras                        = [3.7051, 9.7142, 7.8618, 3.2829, 6.3907, 1.0631, 0.9271, 5.0376, 7.3892, 0.4716, 4.1307, 5.7775, 4.4555, 7.6121, 0.6187, 7.7768, 9.0002, 3.6414, 5.6429, 8.2453, 1.2789, 5.3527];
Model.StepSize                     = 0.05;

% Set up blocking for parameters - either fixed blocks or random blocks
% (This is used for higher-dimensional problems with ill-conditioned metric
% tensors)

% Either 'fixed':
%{
Model.Blocking.Type                = 'fixed';
Model.Blocking.NumOfBlocks         = 1;
Model.Blocking.Blocks{1}           = [1:22];
%Model.Blocking.Blocks{2}           = [2];
%Model.Blocking.Blocks{3}           = [3];
%}

% or 'random'
%
Model.Blocking.Type                = 'random'; % Either 'fixed' or 'random'
Model.Blocking.Blocks              = [4 5 4 5 4]; % Must sum to number of parameters
Model.Blocking.NumOfBlocks         = length(Model.Blocking.Blocks);
%}

% Set up priors for each of the parameters
if Model.Log10Space
    
    % Set up priors for log space
    
    for i = 1:length(Model.Paras)
        Model.Func_Evaluate_Prior{i} = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{i}   = @Sample_Prior_Uniform;
        Model.Priors{i}.Paras        = [-2 2];
    end
    
else
    
    % Set up priors for standard space
    
    for i = 1:length(Model.Paras)
        Model.Func_Evaluate_Prior{i} = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{i}   = @Sample_Prior_Uniform;
        Model.Priors{i}.Paras        = [0 10];
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE model settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up names of ODE solvers for 1st order sensitivities
Model.ModelSpecific.SBModelName_1st_Paras        = 'SBFitzHughNagumo_1st_Paras';

% Set up starting values for all temperatures
% Use log space for initial conditions - helps with very small/large values
Model.ModelSpecific.Log10Space_ICs               = false;
Model.ModelSpecific.ICs                          = [0.1290, 13.6937, 9.1584, 1.9919, 5.9266, 1.1007];

% Set total number of species in ODE system
Model.ModelSpecific.NumOfSpecies                 = 6;

% Set which species are observed and unobserved
Model.ModelSpecific.ObservedSpecies              = 1:6;
Model.ModelSpecific.UnobservedSpecies            = [];

% Set solver tolerances
Model.ModelSpecific.Options.abstol               = 1e-8;
Model.ModelSpecific.Options.reltol               = 1e-8;

% Decide whether to infer initial conditions
Model.ModelSpecific.InferICs                     = true; % otherwise keep them fixed at above values

% Set up blocking for ODE initial conditions - either fixed blocks or random blocks
% (This is used for higher-dimensional problems with ill-conditioned metric
% tensors)
if Model.ModelSpecific.InferICs
    
    % Set up names of ODE solvers for 1st order sensitivities for ICs
    Model.ModelSpecific.SBModelName_1st_ICs              = 'SBFitzHughNagumo_1st_ICs';
    
    % Either 'fixed':
    if strcmp(Model.Blocking.Type, 'fixed')
        
        Model.ModelSpecific.Blocking_ICs.NumOfBlocks     = 2;
        Model.ModelSpecific.Blocking_ICs.Blocks{1}       = [1 2 3];
        Model.ModelSpecific.Blocking_ICs.Blocks{2}       = [4 5 6];
        
        % Must add on the other parameters - these become parameters 4 and 5
        for i = 1:Model.ModelSpecific.Blocking_ICs.NumOfBlocks
            Model.ModelSpecific.Blocking_ICs.Blocks{i} = Model.ModelSpecific.Blocking_ICs.Blocks{i} + Model.NumOfParas;
        end
        
        % Append to current blocking option
        %Model.Blocking.NumOfBlocks = Model.Blocking.NumOfBlocks + Model.ModelSpecific.Blocking_ICs.NumOfBlocks;
        Model.Blocking.Blocks{1}                         = [1:(22+6)];
        
    elseif strcmp(Model.Blocking.Type, 'random')
        % or 'random'
        
        Model.ModelSpecific.Blocking_ICs.Blocks          = [3 3]; % Must sum to number of initial conditions
        Model.ModelSpecific.Blocking_ICs.NumOfBlocks     = length(Model.ModelSpecific.Blocking_ICs.Blocks);
        
        % Append to current blocking option
        Model.Blocking.NumOfBlocks = Model.Blocking.NumOfBlocks + Model.ModelSpecific.Blocking_ICs.NumOfBlocks;
        Model.Blocking.Blocks      = [Model.Blocking.Blocks, Model.ModelSpecific.Blocking_ICs.Blocks];
    end
    %}
    
    % Set up priors for initial conditions
    if Model.ModelSpecific.Log10Space_ICs
        
        % Set up priors for log space
        
        for i = (Model.NumOfParas+1):(Model.NumOfParas+length(Model.ModelSpecific.ICs))
            Model.Func_Evaluate_Prior{i}   = @Evaluate_Prior_Normal;
            Model.Func_Sample_Prior{i}     = @Sample_Prior_Normal;
            Model.Priors{i}.Paras          = [0 1];
        end
        
    else
        
        % Set up priors for standard space
        
        %{
        for i = (Model.NumOfParas+1):(Model.NumOfParas+length(Model.ModelSpecific.ICs))
            %Model.Func_Evaluate_Prior{i}   = @Evaluate_Prior_Uniform;
            %Model.Func_Sample_Prior{i}     = @Sample_Prior_Uniform;
            %Model.Priors{i}.Paras          = [-5 5];
            Model.Func_Evaluate_Prior{i}   = @Evaluate_Prior_Normal;
            Model.Func_Sample_Prior{i}     = @Sample_Prior_Normal;
            Model.Priors{i}.Paras          = [0 10];
        end
        %}
        
        % Prior for IC 1
        Model.Func_Evaluate_Prior{Model.NumOfParas + 1}   = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{Model.NumOfParas + 1}     = @Sample_Prior_Uniform;
        Model.Priors{Model.NumOfParas + 1}.Paras          = [0 0.5];
        
        % Prior for IC 2
        Model.Func_Evaluate_Prior{Model.NumOfParas + 2}   = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{Model.NumOfParas + 2}     = @Sample_Prior_Uniform;
        Model.Priors{Model.NumOfParas + 2}.Paras          = [11 16];
        
        % Prior for IC 3
        Model.Func_Evaluate_Prior{Model.NumOfParas + 3}   = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{Model.NumOfParas + 3}     = @Sample_Prior_Uniform;
        Model.Priors{Model.NumOfParas + 3}.Paras          = [7 12];
       
        % Prior for IC 4
        Model.Func_Evaluate_Prior{Model.NumOfParas + 4}   = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{Model.NumOfParas + 4}     = @Sample_Prior_Uniform;
        Model.Priors{Model.NumOfParas + 4}.Paras          = [0 4];
        
        % Prior for IC 5
        Model.Func_Evaluate_Prior{Model.NumOfParas + 5}   = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{Model.NumOfParas + 5}     = @Sample_Prior_Uniform;
        Model.Priors{Model.NumOfParas + 5}.Paras          = [4 8];
        
        % Prior for IC 6
        Model.Func_Evaluate_Prior{Model.NumOfParas + 6}   = @Evaluate_Prior_Uniform;
        Model.Func_Sample_Prior{Model.NumOfParas + 6}     = @Sample_Prior_Uniform;
        Model.Priors{Model.NumOfParas + 6}.Paras          = [0 3];
        
    end
    
    
    % Append starting initial conditions to parameter vector
    Model.Paras = [Model.Paras Model.ModelSpecific.ICs];
    
    % Infer all initial conditions as well
    Model.NumOfParas = Model.NumOfParas + Model.ModelSpecific.NumOfSpecies;
    
end


%%% Initialise any chain specific values %%%
Model.ChainSpecific = [];


% Set up the data for the ODE model

% Open the given file
InputID                            = 'Locke_Benchmark_Data.mat';
InputData                          = open(InputID);

% Save the data and timepoints
Model.Data                         = InputData.Data;
Model.TimePoints                   = InputData.TimePoints;
NoiseVariance                      = InputData.NoiseVariance;
Model.NumOfTimePoints              = size(Model.TimePoints, 1);

% Save the covariance matrices for calculating the likelihood
for i = Model.ModelSpecific.ObservedSpecies
    Model.NoiseCov{i} = eye(Model.NumOfTimePoints)*NoiseVariance(i);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call main population MCMC routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Population_MCMC(MCMC_Options, Model);



end