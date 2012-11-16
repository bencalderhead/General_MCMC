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
        
        
        
        %%% Set up model specific variables for external functions %%%
        Chain{ChainNum}.ChainSpecific = Model.ChainSpecific;
        
        
        
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
            Chain{ChainNum}.StepSize = ones(1, Model.Blocking.NumOfBlocks)*0.05;
        else
            
            % Use given width
            Chain{ChainNum}.StepSize = ones(1, Model.Blocking.NumOfBlocks)*Model.StepSize;
        end
        
        
        
        %%% Set up proposal counters %%%
        
        Chain{ChainNum}.AcceptedMutation  = zeros(1, Model.Blocking.NumOfBlocks);
        Chain{ChainNum}.AttemptedMutation = zeros(1, Model.Blocking.NumOfBlocks);
        
        Chain{ChainNum}.AcceptedExchange  = 0;
        Chain{ChainNum}.AttemptedExchange = 0;
        
        
        
        %%% Allocate memory for likelihood and geometric quantities %%%
        
        Chain{ChainNum}.LL                    = [];
        Chain{ChainNum}.LogPrior              = [];
        Chain{ChainNum}.GradLL                = [];
        Chain{ChainNum}.GradLogPrior          = [];
        Chain{ChainNum}.FI                    = [];
        Chain{ChainNum}.HessianLogPrior       = [];
        
        
        
        %%% Evaluate the model at new parameters to get LL, gradient, metric etc. %%%
        
        Chain{ChainNum}.CurrentBlock     = [];
        Chain{ChainNum}.CurrentBlockSize = 0;
        Chain{ChainNum}.CurrentBlockNum  = 0;
        
        [Chain{ChainNum}, Success] = Model.Func_Evaluate(Model, Chain{ChainNum});
        
        
        if ~Success
            disp('Error with model evaluation.')
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
        Samples_BurnIn{ChainNum}.LogPrior       = zeros(MCMC_Options.NumOfBurnInSamples, Model.NumOfParas);
        
    end
end



%%% Posterior Samples %%%

for ChainNum = 1:MCMC_Options.NumOfTemps
    
    %%% Allocate memory for posterior samples %%%
    
    Samples_Posterior{ChainNum}.Paras          = zeros(MCMC_Options.PosteriorSaveSize, Model.NumOfParas);
    Samples_Posterior{ChainNum}.LL             = zeros(MCMC_Options.PosteriorSaveSize, 1);
    Samples_Posterior{ChainNum}.LogPrior       = zeros(MCMC_Options.PosteriorSaveSize, Model.NumOfParas);
    
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
    
    %disp(['Iteration no: ' num2str(IterationNum)])
    
    % Repeat for each chain in the population
    parfor ChainNum = 1:MCMC_Options.NumOfTemps
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MCMC update of parameter values
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Check if sampling from the prior
        if Chain{ChainNum}.Temp == 0
            
            
            %%% Sample new values from prior %%%
            
            % Save old chain values
            OldChain = Chain{ChainNum};
            
            % Sample new parameter values
            for n = 1:Model.NumOfParas
                Chain{ChainNum}.Paras(n) = Model.Func_Sample_Prior{n}(Model.Priors{n}.Paras);
            end
            
            % Clear memory for current LL and geometric quantities
            %{
            Chain{ChainNum}.LL                    = [];
            Chain{ChainNum}.LogPrior              = [];
            Chain{ChainNum}.GradLL                = [];
            Chain{ChainNum}.GradLogPrior          = [];
            Chain{ChainNum}.FI                    = [];
            Chain{ChainNum}.HessianLogPrior       = [];
            %}
            
            %%% Evaluate the model at new parameters to get LL, gradient, metric etc. %%%
            Chain{ChainNum}.CurrentBlock     = [];
            Chain{ChainNum}.CurrentBlockSize = 0;
            Chain{ChainNum}.CurrentBlockNum  = 0;
            
            [Chain{ChainNum}, Success] = Model.Func_Evaluate(Model, Chain{ChainNum});
            
            
            if Success
                
            else
                % Reinstate old values if not successful
                Chain{ChainNum} = OldChain;
            end
            
            
        else
            
            if strcmp(Model.Blocking.Type, 'random')
                %%% Set up random blocks %%%
                
                IdxTemp = randperm(Model.NumOfParas);

                for BlockNum = 1:Model.Blocking.NumOfBlocks
                    Start            = sum(Model.Blocking.Blocks(1:BlockNum-1)) + 1;
                    End              = sum(Model.Blocking.Blocks(1:BlockNum));
                    Blocks{BlockNum} = IdxTemp(Start:End);
                end

            else
                %%% Use fixed blocks %%%
                
                Blocks = Model.Blocking.Blocks;
                
            end
            
            
            %%% Update each block of parameters in turn %%%
            
            for BlockNum = 1:Model.Blocking.NumOfBlocks
                
                % Set up variables for updating current block of parameters
                CurrentBlock     = Blocks{BlockNum};
                CurrentBlockSize = length(CurrentBlock);
                
                % Save values in Chain object for passing to functions
                Chain{ChainNum}.CurrentBlock     = CurrentBlock;
                Chain{ChainNum}.CurrentBlockSize = length(CurrentBlock);
                Chain{ChainNum}.CurrentBlockNum  = BlockNum;
                
                
                %%% Sample new parameter values using chosen sampler %%%
                
                % Increment counter
                Chain{ChainNum}.AttemptedMutation(BlockNum) = Chain{ChainNum}.AttemptedMutation(BlockNum) + 1;
                
                
                %%% Evaluate the model at current parameters to get LL, gradient, metric etc. %%%
                
                [Chain{ChainNum}, Success] = Model.Func_Evaluate(Model, Chain{ChainNum});
                
                % Save old chain values
                OldChain = Chain{ChainNum};
                
                
                %%% Calculate the proposal mean and covariance %%%
                try
                    [Proposal_Mean, Proposal_Covariance] = Model.Func_Proposal(Model, Chain{ChainNum});
                catch
                    disp('Error calculating proposal mean and covariance.')
                end
                
                %%% Propose new parameters %%%
                
                RandVec                               = randn(1, CurrentBlockSize);
                Chain{ChainNum}.Paras( CurrentBlock ) = Proposal_Mean + (RandVec*chol(Proposal_Covariance));
                %Chain{ChainNum}.Paras  = mvnrnd(Proposal_Mean, Proposal_Covariance);
                
                % Clear memory for current LL and geometric quantities
                %{
                Chain{ChainNum}.LL                    = [];
                Chain{ChainNum}.LogPrior              = [];
                Chain{ChainNum}.GradLL                = [];
                Chain{ChainNum}.GradLogPrior          = [];
                Chain{ChainNum}.FI                    = [];
                Chain{ChainNum}.HessianLogPrior       = [];
                %}
                
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
                    
                    Prob_NewGivenOld = Normal_LogPDF(Chain{ChainNum}.Paras(CurrentBlock)', Proposal_Mean', Proposal_Covariance);
                    
                    
                    
                    %%% Calculate mean, covariance of proposal and propose new parameters %%%
                    
                    [Proposal_Mean, Proposal_Covariance] = Model.Func_Proposal(Model, Chain{ChainNum});
                    
                    
                    
                    %%% Calculate old given new %%%
                    
                    Prob_OldGivenNew = Normal_LogPDF(OldChain.Paras(CurrentBlock)', Proposal_Mean', Proposal_Covariance);
                    
                    
                    
                    %%% Accept or reject according to ratio %%%
                    
                    Ratio = Chain{ChainNum}.LL*Chain{ChainNum}.Temp + sum(Chain{ChainNum}.LogPrior(CurrentBlock)) + Prob_OldGivenNew - OldChain.LL*Chain{ChainNum}.Temp - sum(OldChain.LogPrior(CurrentBlock)) - Prob_NewGivenOld;
                    
                    
                    
                    if Ratio > 0 || log(rand) < min(0, Ratio)% = log(1) !!
                        % Accept proposal
                        Chain{ChainNum}.AcceptedMutation(BlockNum) = Chain{ChainNum}.AcceptedMutation(BlockNum) + 1;
                        
                    else
                        % Move rejected so use previous chain values
                        Chain{ChainNum} = OldChain;
                        
                        %disp('Rejected')
                    end   
                end
                
                
            end % End block loop
            
        end
        
        
        
        
        
    end % End chain loop
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%% Exchange the chains between temperatures %%%
    if mod(IterationNum, 1) == 0
        
        % Exchange chains...
        for LoopNum = 1:3
            % Loop through each chain in turn and attemp exchange with
            % chain above.
            for ChainNum = 1:(MCMC_Options.NumOfTemps-1)
                
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
                    Chain{Chain2Num}.AttemptedExchange = Chain{Chain1Num}.AttemptedExchange;
                    Chain{Chain2Num}.AcceptedExchange  = Chain{Chain1Num}.AcceptedExchange;
                    Chain{Chain2Num}.StepSize          = Chain{Chain1Num}.StepSize;
                    Chain{Chain2Num}.Temp              = Chain{Chain1Num}.Temp;
                    
                    
                    Chain{Chain1Num}.AttemptedMutation = TempChain.AttemptedMutation;
                    Chain{Chain1Num}.AcceptedMutation  = TempChain.AcceptedMutation;
                    Chain{Chain1Num}.AttemptedExchange = TempChain.AttemptedExchange;
                    Chain{Chain1Num}.AcceptedExchange  = TempChain.AcceptedExchange;
                    Chain{Chain1Num}.StepSize          = TempChain.StepSize;
                    Chain{Chain1Num}.Temp              = TempChain.Temp;
                    
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
            Samples_Posterior{ChainNum}.LogPrior(SavePosition, :)    = Chain{ChainNum}.LogPrior;
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
                Samples_BurnIn{ChainNum}.LogPrior(SavePosition, :)    = Chain{ChainNum}.LogPrior;
            end
            
            
        end
        
        
        
        
        %%% Adjust parameter proposal widths %%%
        
        if mod(IterationNum, MCMC_Options.AdaptRate) == 0
            
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['Burn in iteration: ' num2str(IterationNum) ' of ' num2str(MCMC_Options.NumOfBurnInSamples)]);
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(' ')
            
            % For each population
            for ChainNum = 1:MCMC_Options.NumOfTemps
                
                if (ChainNum > 1) || (MCMC_Options.NumOfTemps == 1)
                
                    for BlockNum = 1:Model.Blocking.NumOfBlocks
                        % Adjust proposal width for parameter value inference

                        AcceptanceRatios(ChainNum, BlockNum) = Chain{ChainNum}.AcceptedMutation(BlockNum)./Chain{ChainNum}.AttemptedMutation(BlockNum);

                        if AcceptanceRatios(ChainNum, BlockNum) < MCMC_Options.Lower
                            Chain{ChainNum}.StepSize(BlockNum) = Chain{ChainNum}.StepSize(BlockNum) * 0.8;
                        elseif AcceptanceRatios(ChainNum, BlockNum) > MCMC_Options.Upper
                            Chain{ChainNum}.StepSize(BlockNum) = Chain{ChainNum}.StepSize(BlockNum) * 1.2;
                        end

                        % Reset counters
                        Chain{ChainNum}.AcceptedMutation(BlockNum)  = 0;
                        Chain{ChainNum}.AttemptedMutation(BlockNum) = 0;
                    end

                end
                
            end
            
            
            %%% Display summary information for each chain %%%
            
            disp('Parameter acceptance rates:')
            
            for BlockNum = 1:Model.Blocking.NumOfBlocks
                for ChainNum = 1:MCMC_Options.NumOfTemps
                    fprintf('%f  ', AcceptanceRatios(ChainNum, BlockNum))
                end
                fprintf('\n')
            end
            fprintf('\n')
            
            disp('Parameter stepsizes:')
            
            for BlockNum = 1:Model.Blocking.NumOfBlocks
                for ChainNum = 1:MCMC_Options.NumOfTemps
                    fprintf('%f  ', Chain{ChainNum}.StepSize(BlockNum))
                end
                fprintf('\n')
            end
            
            fprintf('\n')
            
            
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

save(['./Results/' CurrentFileID], 'Model', 'Chain', 'MCMC_Options', 'Samples_Posterior', 'PosteriorTime');


end


