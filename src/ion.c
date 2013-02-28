#include "ion.h"
#include <fenv.h>

#define MALLOC(m, n, A, lda) \
  do { \
    if ((A = malloc(n * (lda = (m + 1) & ~1) * sizeof(double))) == NULL) \
      return ION_CHANNEL_OUT_OF_MEMORY; \
  } while (0)

#define VALLOC(n, x, incx) \
  do { \
    if ((x = malloc(n * incx * sizeof(double))) == NULL) \
      return ION_CHANNEL_OUT_OF_MEMORY; \
  } while (0)

#define ION_CHANNEL_PROPOSAL_OUTWITH_PRIOR 1

typedef struct {
  double (*pdf)(double, const void *);
  void * params;
} prior;

struct __ion_model_st {
  int n;                // Number of parameters
  prior * priors;       // Prior for each parameter (n of them)
  
  int nStates;     // Number of states
  void (*calculate_Q_matrix)(const double *, double *, size_t);  // Function to generate the Q matrix from the current parameter values
  int nOpenStates, nClosedStates;
  int * openStates, * closedStates;
};

struct __markov_chain_st {
  int n;                // Number of parameters
  double * params;      // Parameters
  double * logprior;    // (log) Prior probability of parameters
}

int ion_evaluate_mh(const ion_model model, markov_chain chain, int * info) {
  // Calculate the prior for each parameter
  for (int i = 0; i < model->n; i++) {
    double x = model->priors[i].pdf(chain->params[i], model->prior[i].params);
    if (x <= 0.0) {
      *info = i + 1;
      return ION_CHANNEL_PROPOSAL_OUTWITH_PRIOR;
    }
    chain->logprior[i] = log(x);
  }


  // Calculate the likelihood of the ion channel data

  // Calculate the Q matrices
  double * Q;
  size_t ldq;
  MALLOC(model->nStates, model->nStates, Q, ldq);
  model->calculate_Q_matrix(chain->params, Q, ldq);

  // Split up the Q matrix into its component matrices  
  double * Q_FF, * Q_FA, * Q_AF, * Q_AA;
  size_t ldqff, ldqfa, ldqaf. ldqaa;
  MALLOC(model->nClosedStates, model->nClosedStates, Q_FF, ldqff);
  MALLOC(model->nClosedStates, model->nOpenStates,   Q_FA, ldqfa);
  MALLOC(model->nOpenStates,   model->nClosedStates, Q_AF, ldqaf);
  MALLOC(model->nOpenStates,   model->nOpenStates,   Q_AA, ldqaa);

  for (int j = 0; j < model->nClosedStates; j++) {
    for (int i = 0; i < model->nClosedStates; i++)
      Q_FF[j * ldqff + i] = Q[model->closedStates[j] * ldq + model->closedStates[i]];
  }
  
  for (int j = 0; j < model->nOpenStates; j++) {
    for (int i = 0; i < model->nClosedStates; i++)
      Q_FA[j * ldqfa + i] = Q[model->openStates[j] * ldq + model->closedStates[i]];
  }
  
  for (int j = 0; j < model->nClosedStates; j++) {
    for (int i = 0; i < model->nOpenStates; i++)
      Q_AF[j * ldqaf + i] = Q[model->closedStates[j] * ldq + model->openStates[i]];
  }
  
  for (int j = 0; j < model->nOpenStates; j++) {
    for (int i = 0; i < model->nOpenStates; i++)
      Q_AA[j * ldqaa + i] = Q[model->openStates[j] * ldq + model->openStates[i]];
  }

  // Calculate equilibrium state occupancies
  // Original code does:
  //  u = ones(1, numstates);
  //  S = [ Q u' ];
  //  EqStates = u / (S * S');
  // This is equivalent to:
  //  u = ones(1, numstates);
  //  eqstates = u / (Q * Q' + ones(numstates));

  double * S, * eqStates;
  size_t lds;
  MALLOC(model->nStates, model->nStates, S, lds);
  VALLOC(model->nStates, eqStates, 1);
  
  // Initialise S to all ones (S = ones(numstates))
  // Initialise eqStates to all ones (eqStates = ones(1, numstates))
  for (int j = 0; j < model->nStates; j++) {
    for (int i = 0; i < model->nStates; i++)
      S[j * lds + i] = 1.0;
    eqStates[j] = 1.0;
  }

  // S += Q * Q' (S = Q * Q' + ones(numstates))
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, model->nStates, model->nStates, model->nStates, 1.0, Q, ldq, Q, ldq, 1.0, S, lds);
  // eqStates /= S (eqStates = ones(1, numstates) / (Q * Q' + ones(numstates)))
  feclearexcept(FE_DIVBYZERO);
  cblas_dgesv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, model->nStates, S, lds, eqStates, 1);      // Numerical instability may occur in here
  if (fegetexcept(FE_DIVBYZERO) != 0)
    return ION_CHANNEL_NUMERICAL_INSTABILITY;
  
  // Don't need Q or S any more
  free(Q);
  free(S);

  double * eqStates_F, * eqStates_A;
  VALLOC(model->nClosedStates, eqStates_F, 1);
  VALLOC(model->nOpenStates, eqStates_A, 1);
  for (int i = 0; i < model->nClosedStates; i++)
    eqStates_F[i] = eqStates[model->closedStates[i]];
  for (int i = 0; i < model->nOpenStates; i++)
    eqStates_A[i] = eqStates[model->openStates[i]];

  // Calculate spectral matrices and eigenvectors of current Q_FF

  // Calculate eigenvectors using QR decomposition in LAPACK
  // If Q_FF is symmetric then X = Y and V_Q_FFi is all zeros (?)
  double * V_Q_FF, * V_QFFi, * X, * Y, * workspace, worksize;
  size_t ldx, ldy;
  VALLOC(model->nClosedStates, V_Q_FF, 1);
  VALLOC(model->nClosedStates, V_Q_FFi, 1);
  MALLOC(model->nClosedStates, model->nClosedStates, X, ldx);
  MALLOC(model->nClosedStates, model->nClosedStates, Y, ldy);
  int lwork = -1, info;
  dgeev_("V", "V", &model->nClosedStates, NULL, &ldq, NULL, NULL, NULL, &ldx, NULL, &ldy, &worksize, &lwork, &info);
  if (info != 0)
    return ION_CHANNEL_EIGENVALUE_FAILED;
  lwork = (int)worksize;
  VALLOC(lwork, workspace, 1);
  dgeev_("V", "V", &model->nClosedStates, Q_FF, &ldqff, V_Q_FF, V_Q_FFi, X, &ldx, Y, &ldy, workspace, &lwork, &info);
  if (info != 0)
    return ION_CHANNEL_EIGENVALUE_FAILED;
  
  free(Q_FF);
  free(V_Q_FFi);
  free(workspace);
  
  // Invert Y using LU decomposition (as it is not symmetric positive definite)
  int * ipiv;
  if ((ipiv = malloc(model->nClosedStates * sizeof(int))) == NULL)
    return ION_CHANNEL_OUT_OF_MEMORY;
  dgetrf_(&model->nClosedStates, &model->nClosedStates, Y, &ldy, ipiv, &info);
  if (info != 0)
    return ION_CHANNEL_INVERSE_FAILED;
  dgetri_(&model->nClosedStates, &model->nClosedStates, Y, &ldy, ipiv, &info);
  if (info != 0)
    return ION_CHANNEL_INVERSE_FAILED;
  
  double ** specMat_Q_FF;
  size_t * ldspecMat_Q_FF;
  if ((specMat_Q_FF = malloc(model->nClosedStates * sizeof(double *))) == NULL)
    return ION_CHANNEL_OUT_OF_MEMORY;
  if ((ldspecMat_Q_FF = malloc(model->nClosedStates * sizeof(size_t))) == NULL)
    return ION_CHANNEL_OUT_OF_MEMORY;
  for (int j = 0; j < model->nClosedStates; j++) {
    MALLOC(model->nClosedStates, model->nClosedStates, specMat_Q_FF[j], ldspecMat_Q_FF[j]);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CBlasTrans, model->nClosedStates, model->nClosedStates, 1, 1.0, &X[j * ldx], 1, &Y[j], ldy, 0.0, specMat_Q_FF[j], ldspecMat_Q_FF[j]);      // Calculate spectral matrices
  }

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
