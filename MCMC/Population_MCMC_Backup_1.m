function [] = Population_MCMC( MCMC_Options, Model )



%%% Start parallel capability if required %%%

matlabpool close force
if MCMC_Options.NumOfCores > 1
    matlabpool('open','local', MCMC_Options.NumOfCores)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise chains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Initialising chains...')
disp(' ')

Chain = cell(1, MCMC_Options.NumOfTemps);


% Initialise each chain in the population
for ChainNum = 1:MCMC_Options.NumOfTemps
    
    % Set success flag in case of problems
    Success = 0;
    
    disp(['Initialising chain ' num2str(ChainNum)])
    
    while ~Success
        
        %%% Set up chain temperature %%%
        Chain{ChainNum}.Temp     = MCMC_Options.Temperatures(ChainNum);
        
        Chain{ChainNum}.Sampler  = Model.Sampler;
        
        
        
        %%% Setup priors and starting parameters %%%
        
        Chain{ChainNum}.Priors = Model.Priors;
        
        if isempty(Model.Paras)
            
            % Randomly sample from prior
            for n = 1:Model.NumOfParas
                Chain{ChainNum}.Paras(n) = Model.Func_Sample_Prior{n}(Model.Priors{n}.Paras);
            end
        else
            
            % Use given parameters
            Chain{ChainNum}.Paras = Model.Paras;
        end
        
        % Set up number of parameters
        Chain{ChainNum}.NumOfParas = Model.NumOfParas;
        
        
        
        %%% Setup starting parameter stepsize %%%
        
        if isempty(Model.StepSize)
            
            % Set parameter step sizes for parameters in each chain in each population
            Chain{ChainNum}.StepSize = 0.05;
        else
            
            % Use given width
            Chain{ChainNum}.StepSize = Model.StepSize;
        end
        
        
        
        %%% Set up proposal counters %%%
        
        Chain{ChainNum}.AcceptedMutation  = 0;
        Chain{ChainNum}.AttemptedMutation = 0;
        
        Chain{ChainNum}.AcceptedExchange  = 0;
        Chain{ChainNum}.AttemptedExchange = 0;
        
        
        %%% Allocate memory for likelihood and geometric quantities %%%
        
        Chain{ChainNum}.LL                    = [];
        Chain{ChainNum}.LogPrior              = [];
        Chain{ChainNum}.GradLL                = [];
        Chain{ChainNum}.GradLogPrior          = [];
        Chain{ChainNum}.FI                    = [];
        Chain{ChainNum}.HessianLogPrior       = [];
        
        
        
        %%% Set up ODE model if required %%%
        
        if strcmp(Model.Type, 'ODE')
            
            if Model.InferICs
                
                %%% Setup priors and starting initial conditions %%%
                
                Chain{ChainNum}.NumOfSpecies = Model.NumOfSpecies;
                
                Chain{ChainNum}.Priors_ICs = Model.Priors_ICs;
                
                if isempty(Model.ICs)
                    
                    % Randomly sample from prior
                    for n = 1:Model.NumOfSpecies
                        Chain{ChainNum}.ICs(n) = Model.Func_Sample_PriorICs{n}(Model.Priors_ICs{n}.Paras);
                    end
                else
                    
                    % Use given intitial conditions
                    Chain{ChainNum}.ICs = Model.ICs;
                end
                
                
                %%% Setup starting intial conditions stepsize %%%
                
                if isempty(Model.StepSize_ICs)
                    
                    % Set parameter step sizes for parameters in each chain in each population
                    Chain{ChainNum}.StepSize_ICs = 0.05;
                else
                    
                    % Use given width
                    Chain{ChainNum}.StepSize_ICs = Model.StepSize_ICs;
                end
                
                
                %%% Set up proposal counters %%%
                
                Chain{ChainNum}.AttemptedMutation_ICs = 0;
                Chain{ChainNum}.AcceptedMutation_ICs  = 0;
                
                
                %%% Allocate memory for likelihood and geometric quantities %%%
                
                Chain{ChainNum}.LogPrior_ICs        = [];
                Chain{ChainNum}.GradLL_ICs          = [];
                Chain{ChainNum}.GradLogPrior_ICs    = [];
                Chain{ChainNum}.FI_ICs              = [];
                Chain{ChainNum}.HessianLogPrior_ICs = [];
                
            else
                
                % Just fix initial conditions
                Chain{ChainNum}.ICs = Model.ICs;
            end
            
        end
        
        
        
        %%% Evaluate the model at current parameters to get LL, gradient, metric etc. %%%
        
        [Chain{ChainNum}, Success] = Model.Func_Evaluate(Model, Chain{ChainNum});
        
        
        
        %%% Evaluate the ODE model at current parameters if required to get LL, gradient, metric etc. %%%
        
        if strcmp(Model.Type, 'ODE')
            if Model.InferICs
                [Chain{ChainNum}, Success] = Model.Func_Evaluate_ICs(Model, Chain{ChainNum});
            end
        end
        
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory to store burn-in and posterior samples as required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('Allocating memory...')
disp(' ')

%%% Burn-in Samples %%%

if MCMC_Options.SaveBurnIn
    
    for ChainNum = 1:MCMC_Options.NumOfTemps
        
        %%% Allocate memory for burn-in samples %%%
        
        Samples_BurnIn{ChainNum}.Paras          = zeros(MCMC_Options.NumOfBurnInSamples, Model.NumOfParas);
        Samples_BurnIn{ChainNum}.LL             = zeros(MCMC_Options.NumOfBurnInSamples, 1);
        Samples_BurnIn{ChainNum}.LogPrior       = zeros(MCMC_Options.NumOfBurnInSamples, 1);
        
        
        %%% Allocate memory for burn-in with ODE model %%%
        
        if strcmp(Model.Type, 'ODE')
            
            if Model.InferICs
                
                Samples_BurnIn{ChainNum}.ICs                = zeros(MCMC_Options.NumOfBurnInSamples, Model.NumOfSpecies);
                Samples_BurnIn{ChainNum}.LogPrior_ICs       = zeros(MCMC_Options.NumOfBurnInSamples, 1);
            else
                
                % If not inferring ICs then just save the current ICs
                Samples_BurnIn{ChainNum}.ICs                = Chain{ChainNum}.ICs;
            end
        end
        
    end
end



%%% Posterior Samples %%%

for ChainNum = 1:MCMC_Options.NumOfTemps
    
    %%% Allocate memory for posterior samples %%%
    
    Samples_Posterior{ChainNum}.Paras          = zeros(MCMC_Options.PosteriorSaveSize, Model.NumOfParas);
    Samples_Posterior{ChainNum}.LL             = zeros(MCMC_Options.PosteriorSaveSize, 1);
    Samples_Posterior{ChainNum}.LogPrior       = zeros(MCMC_Options.PosteriorSaveSize, 1);
    
    
    %%% Allocate memory for posterior with ODE model %%%
    
    if strcmp(Model.Type, 'ODE')
        
        if Model.InferICs
            
            Samples_Posterior{ChainNum}.ICs          = zeros(MCMC_Options.PosteriorSaveSize, Model.NumOfSpecies);
            Samples_Posterior{ChainNum}.LogPrior_ICs = zeros(MCMC_Options.PosteriorSaveSize, 1);
        else
            
            % If not inferring ICs then just save the current ICs
            Samples_Posterior{ChainNum}.ICs          = Chain{ChainNum}.ICs;
        end
    end
    
end



disp('Initialisation Completed.');
disp(' ')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Population MCMC Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Starting population MCMC algorithm...')
disp(' ')



% Main loop
for IterationNum = 1:(MCMC_Options.NumOfBurnInSamples + MCMC_Options.NumOfPosteriorSamples)
    
    
    % Repeat for each chain in the population
    for ChainNum = 1:MCMC_Options.NumOfTemps
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MCMC update of parameter values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Check if sampling from the prior
        if Chain{ChainNum}.Temp == 0
            
            
            %%% Sample new values from prior %%%
            
            % Increment counter
            Chain{ChainNum}.AttemptedMutation = Chain{ChainNum}.AttemptedMutation + 1;
            
            % Save old chain values
            OldChain = Chain{ChainNum};
            
            % Sample new parameter values
            for n = 1:Model.NumOfParas
                Chain{ChainNum}.Paras(n) = Model.Func_Sample_Prior{n}(Model.Priors{n}.Paras);
            end
            
            % Clear memory for current LL and geometric quantities
            Chain{ChainNum}.LL                    = [];
            Chain{ChainNum}.LogPrior              = [];
            Chain{ChainNum}.GradLL                = [];
            Chain{ChainNum}.GradLogPrior          = [];
            Chain{ChainNum}.FI                    = [];
            Chain{ChainNum}.HessianLogPrior       = [];
            
            
            %%% Evaluate the model at new parameters to get LL, gradient, metric etc. %%%
            
            [Chain{ChainNum}, Success] = Model.Func_Evaluate(Model, Chain{ChainNum});
            
            
            if Success
                % Increment counter
                Chain{ChainNum}.AcceptedMutation = Chain{ChainNum}.AcceptedMutation + 1;
            else
                % Reinstate old values if not successful
                Chain{ChainNum} = OldChain;
            end
            
            
        else
            
            
            %%% Sample new parameter values using chosen sampler %%%
            
            % Increment counter
            Chain{ChainNum}.AttemptedMutation = Chain{ChainNum}.AttemptedMutation + 1;
            
            % Save old chain values
            OldChain = Chain{ChainNum};
            
            
            
            %%% Calculate the proposal mean and covariance %%%
            
            [Proposal_Mean, Proposal_Covariance] = Model.Func_Proposal(Chain{ChainNum});
            
            
            
            %%% Propose new parameters %%%
            
            RandVec                = randn(1, Model.NumOfParas);
            Chain{ChainNum}.Paras  = Proposal_Mean + (RandVec*chol(Proposal_Covariance));
            %Chain{ChainNum}.Paras  = mvnrnd(Proposal_Mean, Proposal_Covariance);
            
            % Clear memory for current LL and geometric quantities
            Chain{ChainNum}.LL                    = [];
            Chain{ChainNum}.LogPrior              = [];
            Chain{ChainNum}.GradLL                = [];
            Chain{ChainNum}.GradLogPrior          = [];
            Chain{ChainNum}.FI                    = [];
            Chain{ChainNum}.HessianLogPrior       = [];
            
            
            %%% Evaluate the model at new parameters to get LL, gradient, metric etc. %%%
            
            [Chain{ChainNum}, Success] = Model.Func_Evaluate(Model, Chain{ChainNum});
            
            
            
            %%% Check whether proposed parameter values are valid
            %%% i.e. within prior and likelihood is computable
            if ~Success
                
                % Failed so use previous chain values
                Chain{ChainNum} = OldChain;
                
            else
                
                % Continue
                
                %%% Calculate new given old %%%
                
                Prob_NewGivenOld = Normal_LogPDF(Chain{ChainNum}.Paras', Proposal_Mean', Proposal_Covariance);
                
                
                
                %%% Calculate mean, covariance of proposal and propose new parameters %%%
                
                [Proposal_Mean, Proposal_Covariance] = Model.Func_Proposal(Chain{ChainNum});
                
                
                
                %%% Calculate old given new %%%
                
                Prob_OldGivenNew = Normal_LogPDF(OldChain.Paras', Proposal_Mean', Proposal_Covariance);
                
                
                
                %%% Accept or reject according to ratio %%%
                
                Ratio = Chain{ChainNum}.LL*Chain{ChainNum}.Temp + Chain{ChainNum}.LogPrior + Prob_OldGivenNew - OldChain.LL*Chain{ChainNum}.Temp - OldChain.LogPrior - Prob_NewGivenOld;
                
                
                
                if Ratio > 0 || log(rand) < min(0, Ratio)% = log(1) !!
                    % Accept proposal
                    Chain{ChainNum}.AcceptedMutation = Chain{ChainNum}.AcceptedMutation + 1;
                    
                    if strcmp(Model.Type, 'ODE')
                        if Model.InferICs
                            % Update the gradient, FI etc for the ICs
                            [Chain{ChainNum}, Success] = Model.Func_Evaluate_ICs(Model, Chain{ChainNum});
                        end
                    end
                    
                else
                    % Move rejected so use previous chain values
                    Chain{ChainNum} = OldChain;
                    
                    %disp('Rejected')
                end
                
                
            end
            
            
            
        end
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ODE model - MCMC update of initial conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if strcmp(Model.Type, 'ODE')
            
            if Model.InferICs
                
                
                if Chain{ChainNum}.Temp == 0
                    
                    %%% Sample new values from prior %%%
                    
                    % Increment counter
                    Chain{ChainNum}.AttemptedMutation_ICs = Chain{ChainNum}.AttemptedMutation_ICs + 1;
                    
                    % Save old chain values
                    OldChain = Chain{ChainNum};
                    
                    
                    
                    % Sample new parameter values
                    for n = 1:Model.NumOfSpecies
                        Chain{ChainNum}.ICs(n) = Model.Func_Sample_Prior_ICs{n}(Model.Priors_ICs{n}.Paras);
                    end
                    
                    
                    % Clear memory for current LL and geometric quantities
                    Chain{ChainNum}.LL                        = [];
                    Chain{ChainNum}.LogPrior_ICs              = [];
                    Chain{ChainNum}.GradLL_ICs                = [];
                    Chain{ChainNum}.GradLogPrior_ICs          = [];
                    Chain{ChainNum}.FI_ICs                    = [];
                    Chain{ChainNum}.HessianLogPrior_ICs       = [];
                    
                    
                    %%% Calculate LL, gradient, metric etc. %%%
                    [Chain{ChainNum}, Success] = Model.Func_Evaluate_ICs(Model, Chain{ChainNum});
                    
                    
                    if Success
                        % Increment counter
                        Chain{ChainNum}.AcceptedMutation_ICs = Chain{ChainNum}.AcceptedMutation_ICs + 1;
                    else
                        % Reinstate old values if not successful
                        Chain{ChainNum} = OldChain;
                    end
                    
                    
                else
                    
                    %%% Sample new values using chosen sampler %%%
                    
                    % Increment counter
                    Chain{ChainNum}.AttemptedMutation_ICs = Chain{ChainNum}.AttemptedMutation_ICs + 1;
                    
                    % Save old chain values
                    OldChain = Chain{ChainNum};
                    
                    
                    
                    %%% Calculate the proposal mean and covariance %%%
            
                    [Proposal_Mean, Proposal_Covariance] = Model.Func_Proposal_ICs(Chain{ChainNum});
            
                    
                    
                    %%% Propose new parameters %%%
            
                    RandVec              = randn(1, Model.NumOfSpecies);
                    Chain{ChainNum}.ICs  = Proposal_Mean + RandVec*chol(Proposal_Covariance);

                    
                    % Clear memory for current LL and geometric quantities
                    Chain{ChainNum}.LL                        = [];
                    Chain{ChainNum}.LogPrior_ICs              = [];
                    Chain{ChainNum}.GradLL_ICs                = [];
                    Chain{ChainNum}.GradLogPrior_ICs          = [];
                    Chain{ChainNum}.FI_ICs                    = [];
                    Chain{ChainNum}.HessianLogPrior_ICs       = [];
                    
                    
                    
                    %%% Evaluate the model at new parameters to get LL, gradient, metric etc. %%%
            
                    [Chain{ChainNum}, Success] = Model.Func_Evaluate_ICs(Model, Chain{ChainNum});
            
                    
                    
                    %%% Check whether proposed parameter values are valid
                    %%% i.e. within prior and likelihood is computable
                    if ~Success
                        
                        % Failed so use previous chain values
                        Chain{ChainNum} = OldChain;
                        
                    else
                        
                        % Continue
                        
                        %%% Calculate new given old %%%
                        
                        Prob_NewGivenOld = Normal_LogPDF(Chain{ChainNum}.ICs', Proposal_Mean', Proposal_Covariance);
                        
                        
                        
                        %%% Calculate mean, covariance of proposal and propose new ICs %%%
                
                        [Proposal_Mean, Proposal_Covariance] = Model.Func_Proposal_ICs(Chain{ChainNum});
                
                        
                        
                        %%% Calculate old given new %%%
                
                        Prob_OldGivenNew = Normal_LogPDF(OldChain.ICs', Proposal_Mean', Proposal_Covariance);
                
                        
                        
                        %%% Accept or reject according to ratio %%%
                        
                        Ratio = Chain{ChainNum}.LL*Chain{ChainNum}.Temp + Chain{ChainNum}.LogPrior_ICs + Prob_OldGivenNew - OldChain.LL*Chain{ChainNum}.Temp - OldChain.LogPrior_ICs - Prob_NewGivenOld;
                        
                        
                        if Ratio > 0 || log(rand) < min(0, Ratio)% = log(1) !!
                            % Accept proposal
                            Chain{ChainNum}.AcceptedMutation_ICs = Chain{ChainNum}.AcceptedMutation_ICs + 1;
                            
                            % Update the gradient, FI etc for the parameters
                            [Chain{ChainNum}, Success] = Model.Func_Evaluate(Model, Chain{ChainNum});
                            
                        else
                            % Move rejected so use previous chain values
                            Chain{ChainNum} = OldChain;
                            
                            %disp('Rejected')
                        end
                        
                        
                    end
                    
                    
                    
                end
                
                
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% Exchange the chains between temperatures %%%
    if mod(IterationNum, 1) == 0
        
        % Exchange chains...
        for LoopNum = 1:3
            for ChainNum = 1:(MCMC_Options.NumOfTemps-1)
                
                % Pick two neighbouring chains at random
                %Chain1Num = ceil(rand*(NumOfTemps-1));
                Chain1Num = ChainNum;
                Chain2Num = Chain1Num + 1;
                
                
                % Increment counters
                Chain{Chain1Num}.AttemptedExchange = Chain{Chain1Num}.AttemptedExchange + 1;
                Chain{Chain2Num}.AttemptedExchange = Chain{Chain2Num}.AttemptedExchange + 1;
                
                NewChain1LL = Chain{Chain2Num}.LL;
                NewChain2LL = Chain{Chain1Num}.LL;
                
                NewLL = NewChain1LL*Chain{Chain1Num}.Temp + NewChain2LL*Chain{Chain2Num}.Temp;
                OldLL = Chain{Chain1Num}.LL*Chain{Chain1Num}.Temp + Chain{Chain2Num}.LL*Chain{Chain2Num}.Temp;
                
                Ratio = NewLL - OldLL;
                
                
                if Ratio > 0 || log(rand) < min(0, Ratio)
                    
                    % Update counters
                    Chain{Chain1Num}.AcceptedExchange = Chain{Chain1Num}.AcceptedExchange + 1;
                    Chain{Chain2Num}.AcceptedExchange = Chain{Chain2Num}.AcceptedExchange + 1;
                    
                    % Exchange chains
                    TempChain = Chain{Chain1Num};
                    
                    Chain{Chain1Num} = Chain{Chain2Num};
                    Chain{Chain2Num} = TempChain;
                    
                    % Don't swap these values though!
                    Chain{Chain2Num}.AttemptedMutation = Chain{Chain1Num}.AttemptedMutation;
                    Chain{Chain2Num}.AcceptedMutation  = Chain{Chain1Num}.AcceptedMutation;
                    Chain{Chain2Num}.StepSize          = Chain{Chain1Num}.StepSize;
                    Chain{Chain2Num}.Temp              = Chain{Chain1Num}.Temp;
                    
                    Chain{Chain1Num}.AttemptedMutation = TempChain.AttemptedMutation;
                    Chain{Chain1Num}.AcceptedMutation  = TempChain.AcceptedMutation;
                    Chain{Chain1Num}.StepSize          = TempChain.StepSize;
                    Chain{Chain1Num}.Temp              = TempChain.Temp;
                    
                    if strcmp(Model.Type, 'ODE')
                        if Model.InferICs
                            
                            Chain{Chain2Num}.AttemptedMutation_ICs = Chain{Chain1Num}.AttemptedMutation_ICs;
                            Chain{Chain2Num}.AcceptedMutation_ICs  = Chain{Chain1Num}.AcceptedMutation_ICs;
                            Chain{Chain2Num}.StepSize_ICs          = Chain{Chain1Num}.StepSize_ICs;

                            Chain{Chain1Num}.AttemptedMutation_ICs = TempChain.AttemptedMutation_ICs;
                            Chain{Chain1Num}.AcceptedMutation_ICs  = TempChain.AcceptedMutation_ICs;
                            Chain{Chain1Num}.StepSize_ICs          = TempChain.StepSize_ICs;
                        end
                    end
                                
                    
                    
                    
                end
                
            end
        end
        
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    if IterationNum > MCMC_Options.NumOfBurnInSamples
        
        %%% Save posterior samples %%%
        
        SavePosition = IterationNum - MCMC_Options.NumOfBurnInSamples;
        SavePosition = mod(SavePosition, MCMC_Options.PosteriorSaveSize);
        
        if SavePosition == 0
            SavePosition = MCMC_Options.PosteriorSaveSize;
        end
        
        
        for ChainNum = 1:MCMC_Options.NumOfTemps
            % Update posterior history
            Samples_Posterior{ChainNum}.Paras(SavePosition, :)       = Chain{ChainNum}.Paras;
            Samples_Posterior{ChainNum}.LL(SavePosition)             = Chain{ChainNum}.LL;
            Samples_Posterior{ChainNum}.LogPrior(SavePosition)       = Chain{ChainNum}.LogPrior;
        end
        
        
        
        
        if strcmp(Model.Type, 'ODE')
            %%% Save initial conditions if necessary %%%
            if Model.InferICs
                
                for ChainNum = 1:MCMC_Options.NumOfTemps
                    % Update posterior history
                    Samples_Posterior{ChainNum}.ICs(SavePosition, :)       = Chain{ChainNum}.ICs;
                    Samples_Posterior{ChainNum}.LogPrior_ICs(SavePosition) = Chain{ChainNum}.LogPrior_ICs;
                end
                
                
            end
        end
        
        
        
        %%% Save to file when full %%%
        if SavePosition == MCMC_Options.PosteriorSaveSize && ~(IterationNum == MCMC_Options.PosteriorSaveSize)
            
            CurrentFileID = [MCMC_Options.OutputID '_Posterior_' num2str(IterationNum - MCMC_Options.NumOfBurnInSamples)];
            save(['./Results/' CurrentFileID], 'Model', 'MCMC_Options', 'Samples_Posterior');
            
        end
        
        
        
        %%% Stop if required number of samples has been collected
        if IterationNum == MCMC_Options.NumOfBurnInSamples + 1
            
            % Save burn-in if required
            if MCMC_Options.SaveBurnIn
                CurrentFileID = [MCMC_Options.OutputID '_BurnIn_' num2str(MCMC_Options.NumOfBurnInSamples)];
                save(['./Results/' CurrentFileID], 'Model', 'MCMC_Options', 'Samples_BurnIn');
            end
            
        end
        
        
        if mod(IterationNum, MCMC_Options.AdaptRate) == 0
            disp(['Posterior sample: ' num2str(IterationNum-MCMC_Options.NumOfBurnInSamples) ' of ' num2str(MCMC_Options.NumOfPosteriorSamples)])
        end
        
        
        
    else
        
        if MCMC_Options.SaveBurnIn
            
            %%% Save burn-in samples if required %%%
            
            SavePosition = IterationNum;
            
            
            for ChainNum = 1:MCMC_Options.NumOfTemps
                % Update posterior history
                Samples_BurnIn{ChainNum}.Paras(SavePosition, :)       = Chain{ChainNum}.Paras;
                Samples_BurnIn{ChainNum}.LL(SavePosition)             = Chain{ChainNum}.LL;
                Samples_BurnIn{ChainNum}.LogPrior(SavePosition)       = Chain{ChainNum}.LogPrior;
            end
            
            
            if strcmp(Model.Type, 'ODE')
                %%% Save initial conditions if necessary %%%
                if Model.InferICs
                    
                    for ChainNum = 1:MCMC_Options.NumOfTemps
                        % Update posterior history
                        Samples_BurnIn{ChainNum}.ICs(SavePosition, :)       = Chain{ChainNum}.ICs;
                        Samples_BurnIn{ChainNum}.LogPrior_ICs(SavePosition) = Chain{ChainNum}.LogPrior_ICs;
                    end
                    
                end
            end
            
            
        end
        
        
        
        
        %%% Adjust parameter proposal widths %%%
        
        if mod(IterationNum, MCMC_Options.AdaptRate) == 0
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['Burn in iteration: ' num2str(IterationNum) ' of ' num2str(MCMC_Options.NumOfBurnInSamples)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(' ')
            
            % For each population
            for ChainNum = 2:MCMC_Options.NumOfTemps
                
                % Adjust proposal width for parameter value inference
                
                AcceptanceRatios(ChainNum) = Chain{ChainNum}.AcceptedMutation./Chain{ChainNum}.AttemptedMutation;
                
                if AcceptanceRatios(ChainNum) < MCMC_Options.Lower
                    Chain{ChainNum}.StepSize = Chain{ChainNum}.StepSize * 0.8;
                elseif AcceptanceRatios(ChainNum) > MCMC_Options.Upper
                    Chain{ChainNum}.StepSize = Chain{ChainNum}.StepSize * 1.2;
                end
                
                % Reset counters
                Chain{ChainNum}.AcceptedMutation  = 0;
                Chain{ChainNum}.AttemptedMutation = 0;
                
                
                if strcmp(Model.Type, 'ODE')
                    
                    % Adjust proposal width for inital condition inference
                    
                    if Model.InferICs
                        
                        AcceptanceRatiosICs(ChainNum) = Chain{ChainNum}.AcceptedMutation_ICs./Chain{ChainNum}.AttemptedMutation_ICs;
                        
                        if AcceptanceRatiosICs(ChainNum) < Model.Lower_ICs
                            Chain{ChainNum}.StepSize_ICs = Chain{ChainNum}.StepSize_ICs * 0.8;
                        elseif AcceptanceRatiosICs(ChainNum) > Model.Upper_ICs
                            Chain{ChainNum}.StepSize_ICs = Chain{ChainNum}.StepSize_ICs * 1.2;
                        end
                        
                        % Reset counters
                        Chain{ChainNum}.AcceptedMutation_ICs  = 0;
                        Chain{ChainNum}.AttemptedMutation_ICs = 0;
                        
                    end
                end
                
            end
            
            
            %%% Display summary information for each chain %%%
            
            disp('Parameter acceptance rates:')
            for ChainNum = 1:MCMC_Options.NumOfTemps
                fprintf('%f  ', AcceptanceRatios(ChainNum))
            end
            fprintf('\n')
            
            disp('Parameter stepsizes:')
            for ChainNum = 1:MCMC_Options.NumOfTemps
                fprintf('%f  ', Chain{ChainNum}.StepSize)
            end
            fprintf('\n')
            
            if strcmp(Model.Type, 'ODE')
                if Model.InferICs
                    
                    disp(' ')
                    disp('IC acceptance rates:')
                    for ChainNum = 1:MCMC_Options.NumOfTemps
                        fprintf('%f  ', AcceptanceRatiosICs(ChainNum))
                    end
                    fprintf('\n')
                    
                    disp('IC stepsizes:')
                    for ChainNum = 1:MCMC_Options.NumOfTemps
                        fprintf('%f  ', Chain{ChainNum}.StepSize_ICs)
                    end
                    fprintf('\n')
                    
                end
            end
            
            % Display exchange rate
            disp(' ')
            disp('Model parameter exchange ratios:')
            for ChainNum = 1:MCMC_Options.NumOfTemps
                fprintf('%f  ', Chain{ChainNum}.AcceptedExchange./Chain{ChainNum}.AttemptedExchange)
                
                Chain{ChainNum}.AcceptedExchange  = 0;
                Chain{ChainNum}.AttemptedExchange = 0;
            end
            fprintf('\n')
            fprintf('\n')
            
            
        end
        
        
        
        
        
        
        % Check for simple convergence
        if IterationNum == MCMC_Options.NumOfBurnInSamples
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp(['Burn in procedure completed.']);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
            disp(' ')
            
            % Begin timing
            tic;
        end
        
    end
    
    
end


PosteriorTime = toc;

% Close parallel capability
if MCMC_Options.NumOfCores > 1
    matlabpool close
end

CurrentFileID = [MCMC_Options.OutputID '_Posterior_' num2str(MCMC_Options.NumOfPosteriorSamples) 'Finish'];

save(['./Results/' CurrentFileID], 'Model', 'MCMC_Options', 'Samples_Posterior', 'PosteriorTime');


end


