#include "ion.h"

#define ION_CHANNEL_PROPOSAL_OUTWITH_PRIOR 1

typedef struct {
  double (*pdf)(double, const void *);
  void * params;
} prior;

struct __ion_model_st {
  int n;                // Number of parameters
  prior * priors;       // Priors (n of them)
  double * Q, size_t ldq;           // Q matrix (n by n)
  void (*calculate_Q_matrix)(const double *, double *);  // Function to generate the Q matrix
  int closed, open;           // Offsets of closed and open states to split the Q matrix
};

struct __markov_chain_st {
  unsigned int n;       // Number of parameters
  double * params;      // Parameters
  double * logprior;    // (log) Prior probability of parameters
}

int ion_evaluate_mh(const ion_model model, markov_chain chain, int * info) {
  // Calculate the prior for each parameters
  for (int i = 0; i < model->n; i++) {
    double x = model->priors[i].pdf(chain->params[i], model->prior[i].params);
    if (x > 0.0)
      chain->logprior[i] = log(x);
    else {
      *info = i + 1;
      return ION_CHANNEL_PROPOSAL_OUTWITH_PRIOR;
    }
  }


  // Calculate the likelihood of the ion channel data

  // Calculate the Q matrices
  model->calculate_Q_matrix(chain->params, model->Q, model->ldq);

  // Split up the Q matrix into its component matrices
  double * Q_FF = &model->Q[model->closed * model->ldq + model->closed];
  double * Q_FA = &model->Q[model->closed * model->ldq + model->open];
  double * Q_AF = &model->Q[model->open * model->ldq + model->closed];
  double * Q_AA = &model->Q[model->open * model->ldq + model->open];

  // Calculate equilibrium state occupancies
  // Original code does:
  //  u = ones(1, numstates);
  //  S = [ Q u' ];
  //  EqStates = u / (S * S');
  // This is equivalent to:
  //  eqstates = 1 / (Q * Q' + 1);
  // Not sure if it is possible to do everything in place to avoid (slow) malloc
  // to store Q * Q'
  double * S, * eq_states;
  if ((S = malloc(model->nstates * model->nstates * sizeof(double))) == NULL)
    return ION_CHANNEL_OUT_OF_MEMORY;
  if ((eq_states = malloc(model->nstates * sizeof(double))) == NULL)
    return ION_CHANNEL_OUT_OF_MEMORY;

  for (int j = 0; j < model->nstates; j++) {
    for (int i = 0; i < model->nstates; i++)
      S[j * model->nstates + i] = 1.0;
    eq_states[j] = 1.0;
  }

  dgemm(CblasNoTrans, CblasTrans, model->nstates, model->nstates, model->nstates, 1.0, Q, ldq, Q, ldq, 1.0, S, model->nstates);
  dgesv(CblasUpper, CblasNoTrans, CblasNonUnit, model->nstates, S, model->nstates, eq_states, 1);      // Numerical instability

  free(S);

  double * eq_states_f = &eq_states[model->closed];
  double * eq_states_a = &eq_states[model->open];

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




}

int ion_proposal_mh(const ion_model, const markov_chain, double *, double *);
