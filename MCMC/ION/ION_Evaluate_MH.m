function [ Chain, Success ] = ION_Evaluate_MH( Model, Chain )

% Inputs: Model object, Chain object
% Output: Chain object, Success flag


% Initialise success flag
Success = 1;


NumOfParas   = Model.NumOfParas; % Number of parameters



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



%%% Calculate the likelihood of the ion channel data %%%

% Calculate the Q matrices
if Model.Log10Space
    Q = Calculate_Q_Matrix(Model.Name, 10.^(Chain.Paras)); % Convert from log space
else
    Q = Calculate_Q_Matrix(Model.Name, Chain.Paras);
end


% Split up the Q matrix into its component matrices
Q_FF = Q(Model.ModelSpecific.ClosedStates, Model.ModelSpecific.ClosedStates);
Q_FA = Q(Model.ModelSpecific.ClosedStates, Model.ModelSpecific.OpenStates);
Q_AA = Q(Model.ModelSpecific.OpenStates, Model.ModelSpecific.OpenStates);
Q_AF = Q(Model.ModelSpecific.OpenStates, Model.ModelSpecific.ClosedStates);


% Calculate equilibrium state occupancies

u = ones(1,Model.ModelSpecific.NumOfStates);
S = [Q u'];

lastwarn('')
EqStates = u/(S*S');

[msgstr, msgid] = lastwarn;
if ~isempty(msgstr)
    disp('Possible numerical instability')
    Success = 0;
    return
end


EqStates_F = EqStates(Model.ModelSpecific.ClosedStates);
EqStates_A = EqStates(Model.ModelSpecific.OpenStates);



% Calculate spectral matrices and eigenvectors of current Q_FF
[X V] = eig(Q_FF);
V_Q_FF = diag(V); % Eigenvectors
Y = inv(X);
for j = 1:length(EqStates_F)
    SpecMat_Q_FF{j} = X(:,j)*Y(j,:); % Calculate spectral matrices
end

% Calculate spectral matrices and eigenvectors of current Q_AA
[X V] = eig(Q_AA);
V_Q_AA = diag(V); % Eigenvectors
Y = inv(X);
for j = 1:length(EqStates_A)
    SpecMat_Q_AA{j} = X(:,j)*Y(j,:); % Calculate spectral matrices
end


% Calculate initial vectors for current state
if Model.Data(1) == 0 % if starts with closed state
    
    L = EqStates_F; % Equilibrium closed states
else
    
    L = EqStates_A; % Equilibrium open states
end

% Calculate in log-space
LL = log(L);


% Do in a slow loop to begin with
for i = 1:(Model.NumOfTimePoints-1)
    
    if Model.Data(i) == 0 % Currently in closed state (denoted F in literature)
        
        % Moving to open state
        Sojourn = Model.TimePoints(i+1) - Model.TimePoints(i); % Get time interval to next move
        
        % Calculate idealised transition probability from closed to open
        
        G_FA = zeros(length(EqStates_F));
        
        for j = 1:length(EqStates_F)
            G_FA = G_FA + exp( V_Q_FF(j)*Sojourn ).*SpecMat_Q_FF{j};
        end
        
        G_FA = G_FA*Q_FA;
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Logarithmic update %
        %%%%%%%%%%%%%%%%%%%%%%
        % Calculate log(L*G_FA) in terms of LL = log(L) and log(G_FA)
        Log_G_FA = log(G_FA);
        
        % Do element-wise logarithmic calculation
        % log(vector*matrix) in terms of log(vector) and log(matrix)
        %
        LL = dlngevm(LL, Log_G_FA);
        %LL(isnan(LL)) = -inf;
        %}
        
        %{
        NewLL = zeros(1, length(EqStates_A)); % Reset to zeros
        
        for j = 1:length(EqStates_A)
            % Multiply row by column in log - i.e. sum!
            S = LL + Log_G_FA(:, j)';
            
            if length(S) > 1
                % Now add them all together in log to get jth element
                S = sort(S, 'descend'); % Sort with biggest first
                
                % Check the first values is not -inf
                if S(1) ~= -inf % Then do normal exponential sum
                    NewLL(j) = S(1) + log( 1 + sum( exp(S(2:end) - S(1)) ) );
                else
                    NewLL(j) = -inf;
                end
            else
                % Only one value so no need to sum!
                NewLL(j) = S;
            end
        end
        
        % Copy new LL
        LL = NewLL;
        %}
        
    else % Currently in open state (denoted A in literature)
        
        % Moving to closed state
        Sojourn = Model.TimePoints(i+1) - Model.TimePoints(i); % Get time interval to next move
        
        % Calculate idealised transition probability from closed to open
        
        G_AF = zeros(length(EqStates_A));
        
        for j = 1:length(EqStates_A)
            G_AF = G_AF + exp( V_Q_AA(j)*Sojourn ).*SpecMat_Q_AA{j};
        end
        
        G_AF = G_AF*Q_AF;
        
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Logarithmic update %
        %%%%%%%%%%%%%%%%%%%%%%
        % Calculate log(L*G_FA) in terms of LL = log(L) and log(G_FA)
        Log_G_AF = log(G_AF);
        
        % Do element-wise logarithmic calculation
        % log(vector*matrix) in terms of log(vector) and log(matrix)
        
        %
        LL = dlngevm(LL, Log_G_AF);
        %LL(isnan(LL)) = -inf;
        %}
        
        %{
        NewLL = zeros(1, length(EqStates_F)); % Reset to zeros
        
        for j = 1:length(EqStates_F)
            % Multiply row by column in log - i.e. sum!
            S = LL + Log_G_AF(:, j)';
            
            if length(S) > 1
                % Now add them all together in log to get jth element
                S = sort(S, 'descend'); % Sort with biggest first
                
                % Check the first values is not -inf
                if S(1) ~= -inf % Then do normal exponential sum
                    NewLL(j) = S(1) + log( 1 + sum( exp(S(2:end) - S(1)) ) );
                else
                    NewLL(j) = -inf;
                end
            else
                % Only one value so no need to sum!
                NewLL(j) = S;
            end
        end
        
        % Copy new LL
        LL = NewLL;
        %}
    end
    
    
end



% Actually all this unit vector multiplication does is sum the likelihood

% Sum the log-likelihood terms
if length(LL) > 1
    % Now add them all together in log to get jth element
    S = sort(LL, 'descend'); % Sort with biggest first
    
    % Check first value is not -inf!
    if S(1) > -inf
        LL = S(1) + log( 1 + sum( exp(S(2:end) - S(1)) ) );
    else
        LL = -inf;
    end
else
    % Only one value so no need to sum!
end


Chain.LL = LL;


if isnan(Chain.LL) || isinf(Chain.LL)
    Success = 0;
    Chain.LL = -1e300;
end




end

