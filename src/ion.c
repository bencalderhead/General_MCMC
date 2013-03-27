#include <gmcmc/ion.h>#ifndef GMCMC_ION_H
#define GMCMC_ION_H

#include <gmcmc/error.h>

/*!
 * ION model structure.
 */
typedef struct __ion_model_st {
  int open, closed;                                     /*!< Number of open and closed states */
  void (*update_q)(const double *, double *, size_t);   /*!< Function to update Q matrix based on current parameter values */
} gmcmc_ion_model;

/**
 * ION model type.
 */
extern const gmcmc_model_type * gmcmc_ion_model_type;

gmcmc_error gmcmc_ion_model_create(gmcmc_ion_model **, int, int, void (*)(const double *, double *, size_t));
void gmcmc_ion_model_destroy(gmcmc_ion_model *);

/*!
 * Calculates the log likelihood of the ion channel data for the current parameter values.
 *
 * @param [in]  model  the model
 * @param [out] chain  the chain
 *
 * @return GMCMC_SUCCESS on success.
 */
gmcmc_error gmcmc_ion_evaluate_mh(const gmcmc_model *, gmcmc_chain *);

/**
 * Calculates the proposal mean and covariance based on the current parameter
 * values.
 *
 * @param [in]  model  the model
 * @param [in]  chain  the chain
 * @param [out] mean   the proposal mean
 * @param [out] cov    the proposal covariance
 *
 * @return GMCMC_SUCCESS on success.
 */
gmcmc_error gmcmc_ion_proposal_mh(const gmcmc_model *, const gmcmc_chain *, double *, double *);

#endif


#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fenv.h>

#include "mcmc_error.h"

/**
 * Calculates the eigenvectors and eigenvalues of an n-by-n real nonsymmetric
 * matrix A.
 *
 * The computed eigenvectors are normalized to have Euclidean norm equal to 1
 * and largest component real.
 *
 * @param n    the size of the matrix A.
 * @param A    a pointer to an array of size n * ldx containing the elements of
 *               A in column major layout.  A is overwritten on output.
 * @param lda  the leading dimension of the array A.
 * @param B    a pointer to an array of size n * ldb.  The eigenvectors of A are
 *               stored column by column in B.
 * @param ldb  the leading dimension of the array B.
 * @param x    a pointer to an array of size n.  The eigenvalues of A are stored
 *               in x.
 */
static int eigen(int, double *, int, double *, int, complex *);
static int inv(int, const double *, int, double *, int);

struct __ion_workspace_st {
  // Q matrix
  double * Q;
  size_t ldq;
  // Equilibrium state occupancies
  double * eqstates;
  // Eigenvectors of Q_FF and Q_AA
  double * v_Q_FF, * v_Q_AA;
  // Spectral matrices of Q_FF and Q_AA
  double * specMat_Q_FF, * specMat_Q_AA;
  size_t ldqff, ldqaa;
  // Temporary matrices (nstates by nstates)
  double * X, * Y;
  // Workspace for the eigenvector and inverse functions
  double * work;
  int * ipiv, lwork;
}

static void ion_workspace_destroy(ion_workspace * workspace) {
  free(workspace->Q);
  free(workspace->eqstates);
  free(workspace->v_Q_FF);
  free(workspace->v_Q_AA);
  free(workspace->specMat_Q_FF);
  free(workspace->specMat_Q_AA);
  free(workspace->X);
  free(workspace->Y);
  free(workspace->work);
  free(workspace->ipiv);
  free(workspace);
}

static mcmc_error ion_workspace_create(ion_workspace ** workspace, unsigned int open, unsigned int closed) {
  if ((*workspace = calloc(1, sizeof(ion_workspace))) == NULL)
    return MCMC_ERROR_OUT_OF_MEMORY;

  // Total number of states
  unsigned int n = open + closed;

  // Allocate the Q matrix
  (*workspace)->ldq = (n + 1u) & ~1u;
  if (((*workspace)->Q = malloc(n * (*workspace)->ldq * sizeof(double))) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate the equilibrium state vector
  if (((*workspace)->eqstates = malloc(n * sizeof(double))) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate vectors for the eigenvalues of Q_FF and Q_AA
  if (((*workspace)->v_Q_FF = malloc(closed * sizeof(double))) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }
  if (((*workspace)->v_Q_AA = malloc(open * sizeof(double))) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate the spectral matrices
  (*workspace)->ldqff = (closed + 1u) & ~1u;
  (*workspace)->ldqaa = (open + 1u) & ~1u;
  if (((*workspace)->specMat_Q_FF = malloc(closed * closed * (*workspace)->ldqff)) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }
  if (((*workspace)->specMat_Q_AA = malloc(open * open * (*workspace)->ldqaa)) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate temporary matrices
  (*workspace)->ldx = (n + 1u) & ~1u;
  (*workspace)->ldy = (n + 1u) & ~1u;
  if (((*workspace)->X = malloc(n * (*workspace)->ldx * sizeof(double))) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }
  if (((*workspace)->Y = malloc(n * (*workspace)->ldy * sizeof(double))) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  // Work out the maximum workspace needed to calculate the eigenvectors and inverse
  int lwork = -1, info;
  double eigensize;
  dgeev_("N", "V", &n, NULL, &n, NULL, NULL, NULL, &n, NULL, &n, &eigensize, &lwork, &info);
  if (info != 0) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_INVALID_ARGUMENT;
  }
  double invsize;
  dgetri_(&n, NULL, &ldy, NULL, &invsize, &lwork, &info);
  if (info != 0) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_INVALID_ARGUMENT;
  }
  if (invsize > eigensize) {
    (*workspace)->lwork = (int)invsize;
    (*workspace)->workspace = malloc((size_t)invsize * sizeof(double));
  }
  else {
    (*workspace)->lwork = (int)eigensize;
    (*workspace)->workspace = malloc((size_t)eigensize * sizeof(double));
  }
  if (((*workspace)->ipiv = malloc(n * sizeof(int))) == NULL) {
    ion_workspace_destroy(*workspace);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  return MCMC_SUCCESS;
}

int ion_model_create(ion_model * restrict model, const char * name, unsigned int open, unsigned int closed, void (*update_q)(const double *, double *, size_t)) {
  mcmc_error error = mcmc_model_create(&model->m);
  if (error != 0)
    return error;

  model->open = open;
  model->closed = closed;
  model->update_q = update_q;

  return ion_workspace_create(&model->workspace, open, closed);
}

void ion_model_destroy(ion_model * model) {
  ion_workspace_destroy(model->workspace);
  mcmc_model_destroy(&model->m);
}

mcmc_error ion_evaluate_mh(const ion_model * restrict model, markov_chain * restrict chain) {
  // Calculate the prior for each parameter
  for (int i = 0; i < model->n; i++) {
    double x = prior_evaluate(&model->priors[i], chain->params[i]);
    if (isnan(x) || islessequal(x, 0.0))
      return MCMC_ERROR_PROPOSAL_OUTWITH_PRIOR;
    chain->logprior[i] = log(x);
  }

  // Calculate the likelihood of the ion channel data
  const unsigned int nstates = model->open + model->closed;

  // Calculate the Q matrix
  size_t ldq = model->workspace->ldq;
  double * Q = model->workspace->Q;
  model->update_q(chain->params, Q, ldq);

  // Calculate equilibrium state occupancies
  double * eqstates = model->workspace->eqstates;
  MCMC_ERROR_CHECK(calculate_equilibrium_state_occupancies(nstates, Q, ldq, eqstates, model->workspace));

  // Split up the Q matrix into its component matrices
  double * Q_FF = Q;
  double * Q_FA = &Q[model->open * ldq];
  double * Q_AF = &Q[model->open];
  double * Q_AA = &Q[model->open * ldq + model->open];

  // Split up the equilibrium states
  double * eqstates_F = eqstates;
  double * eqstates_A = &eqstates[model->open];


  // Calculate spectral matrices and eigenvectors of current Q_FF
  size_t ldqff = model->workspace->ldqff;
  double * v_Q_FF = model->workspace->v_Q_FF;
  double * specMat_Q_FF = model->workspace->specMat_Q_FF;
  MCMC_ERROR_CHECK(calculate_specmat_eigenvectors(model->closed, Q_FF, ldq, specMat_Q_FF, ldqff, v_Q_FF, model->workspace));

  // Calculate spectral matrices and eigenvectors of current Q_AA
  size_t ldqaa = model->workspace->ldqaa;
  double * v_Q_AA = model->workspace->v_Q_AA;
  double * specMat_Q_FF = model->workspace->specMat_Q_AA;
  MCMC_ERROR_CHECK(calculate_specmat_eigenvectors(model->closed, Q_AA, ldq, specMat_Q_AA, ldqaa, v_Q_AA, model->workspace));


  // Calculate initial vectors for current state
  double * L =  (model->data[0] == 0) ? eqStates_F : eqStates_A;


  mcmc_model * m = model->m;
  // Do in a slow loop to begin with
  for (int i = 0; i < m->ntimepoints - 1; i++) {
    double sojourn = m->timepoints[i + 1] - m->timepoints[i]; // Get time interval to next move

    if (m->data[i] == 0) {   // Currently in closed state (denoted F in literature)

      // Calculate idealised transition probability from closed to open

      size_t ldgfa = model->workspace->ldx;
      double * G_FA = model->workspace->X;

      for (int j = 0; j < model->closed; j++) {
        for (int i = 0; i < model->closed; i++)
          G_FA[j * ldgfa + i] = 0.0;
      }

      for (int k = 0; k < model->closed; k++) {
        double alpha = exp(V_Q_FF[j] * sojourn);
        double * specMat = &specMat_Q_FF[k * model->closed * ldsqff];
        for (int j = 0; j < model->closed; j++) {
          for (int i = 0; i < model->closed; i++)
            G_FA[j * lda + i] += alpha * specMat[j * ldsqff + i];
        }
      }

      double * Y = model->workspace->Y;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, model->closed, model->closed, model->closed, 1.0, G_FA, ldgfa, Q_FA, ldqfa, 0.0, Y, ldgfa);
      G_FA = Y;

      // Logarithmic update
      // Calculate log(L*G_FA) in terms of LL = log(L) and log(G_FA)
      // Do element-wise logarithmic calculation
      // log (vector * matrix) in terms of log(vector) and log(matrix)
      log_dgemv(CblasTrans, 1.0, G_FA, ldgfa, L, 1, 0.0, L, 1);
    }
    else {      // Currently in open state (denoted A in literature)

      // Moving to open state
      double sojourn = model->timePoints[i + 1] - model->timePoints[i]; // Get time interval to next move

      // Calculate idealised transition probability from closed to open

      double * G_AF;
      int ldgaf;
      ERROR_CHECK(matrix_alloc(model->nOpenStates, model->nOpenStates, (void **)&G_AF, &ldgaf, sizeof(double)));
      matrix_init(model->nOpenStates, model->nOpenStates, 0.0, G_AF, ldgaf);

      for (int j = 0; j < model->nOpenStates; j++) {
        double alpha = exp(V_Q_AA[j] * sojourn);
        for (int k = 0; k < model->nOpenStates; k++)
          cblas_daxpy(model->nOpenStates, alpha, &SpecMat_Q_AA[j][k * ldspecMat_Q_AA[j]], 1, &G_AF[k * ldgaf], 1);
      }

      double * T;
      int ldt;
      ERROR_CHECK(matrix_alloc(model->nOpenStates, model->nOpenStates, (void **)&T, &ldt, sizeof(double)));
      matrix_copy(model->nOpenStates, model->nOpenStates, G_FA, ldgfa, T, ldt);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, model->nOpenStates, model->nOpenStates, model->nOpenStates, 1.0, T, ldt, Q_AF, ldqfa, 1.0, G_AF, ldqaf);
      free(T);

      // Logarithmic update
      // Calculate log(L*G_FA) in terms of LL = log(L) and log(G_FA)
      Log_G_AF = log(G_AF);

      // Do element-wise logarithmic calculation
      // log (vector * matrix) in terms of log(vector) and log(matrix)
      LL = dlngevm(LL, Log_G_AF);

    }
  }

  // Actually all this unit vector multiplication does is sum the likelihood

  // Sum the log-likelihood terms
  if (length(LL) > 1) {
      // Now add them all together in log to get jth element
      S = sort(LL, 'descend'); % Sort with biggest first

      % Check first value is not -inf!
      if S(1) > -inf
          LL = S(1) + log( 1 + sum( exp(S(2:end) - S(1)) ) );
      else
          LL = -inf;
      end
  }


Chain.LL = LL;


if isnan(Chain.LL) || isinf(Chain.LL)
    Success = 0;
    Chain.LL = -1e300;
end




}

int ion_proposal_mh(const ion_model, const markov_chain, double *, double *);

static int eig(int n, const double * X, int ldx, double complex * V, int ldv, double complex * D, int ldd) {
  if (n == 0)
    return 0;

  // Calling DGEEV with lwork = -1 causes it to return the optimal workspace size in work[0]
  double size;
  int lwork = -1, info;
  dgeev_("N", "V", &n, NULL, &n, NULL, NULL, NULL, &n, NULL, &n, &size, &lwork, &info);
  if (info != 0)
    return info;

  // The real and imaginary components of the eigenvalues are stored separately in wr and wi
  // The (right) eigenvectors are stored column by column in VR
  double * A, * wr, * wi, * VR, * work;
  int lda, ldvl = 1, ldvr;   // ldvl must be >= 1 even when jobvl = "N"
  lwork = (int)size;
  ERROR_CHECK(matrix_alloc(n, n, (void **)&A, &lda, sizeof(double)));
  ERROR_CHECK(vector_alloc(n, (void **)&wr, 1, sizeof(double)));
  ERROR_CHECK(vector_alloc(n, (void **)&wi, 1, sizeof(double)));
  ERROR_CHECK(matrix_alloc(n, n, (void **)&VR, &lda, sizeof(double)));
  ERROR_CHECK(vector_alloc(lwork, (void **)&work, 1, sizeof(double)));

  // Copy X into A as X would be overwritten by the in-place LAPACK routine
  matrix_copy(n, n, NULL, NULL, X, ldx, A, lda);

  // Calculate the eigenvalues and (right) eigenvectors
  dgeev_("N", "V", &n, A, &lda, wr, wi, NULL, &ldvl, VR, &ldvr, work, &lwork, &info);

  // V is an n by n matrix with an eigenvector in each column
  for (int j = 0; j < n; j++) {
    // If the j-th eigenvalue is real, then v(j) = VR(:,j), the j-th column of VR.
    if (wi[j] == 0.0) {
      for (int i = 0; i < n; i++)
        V[j * ldv + i] = VR[j * ldvr + i];
    }
    // If the j-th and (j+1)-st eigenvalues form a complex conjugate pair, then
    // v(j) = VR(:,j) + i*VR(:,j+1) and v(j+1) = VR(:,j) - i*VR(:,j+1).
    else if (wi[j] == -wi[j + 1]) {
      double * a = &VR[j * ldvr], * b = &VR[(j + 1) * ldvr];
      for (int i = 0; i < n; i++)
        V[j * ldv + i] = a[i] + b[i] * I;
      j++;
      for (int i = 0; i < n; i++)
        V[j * ldv + i] = a[i] - b[i] * I;
    }
  }

  // D is a diagonal matrix with the eigenvalues along the diagonal
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < j; i++)
      D[j * ldd + i] = 0.0;
    D[j * ldd + j] = wr[j] + wi[j] * I;
    for (int i = j + 1; i < n; i++)
      D[j * ldd + i] = 0.0;
  }

  free(A);
  free(wr);
  free(wi);
  free(VR);
  free(work);

  return info;
}

static int inv(int n, const double * X, int ldx, double * Y, int ldy) {
  if (n == 0)
    return 0;

  // Copy X into Y as LAPACK is in place
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++)
      Y[j * ldy + i] = X[j * ldx + i];
  }

  // Pivot array
  int * ipiv;
  if ((ipiv = malloc((size_t)n * sizeof(int))) == NULL)
    return ENOMEM;

  // Calculate the LU decomposition of Y
  int info;
  dgetrf_(&n, &n, Y, &ldy, ipiv, &info);

  if (info != 0)
    return info;

  // Calling DGETRI with lwork = -1 causes it to return the optimal workspace size in work[0]
  double size;
  int lwork = -1;
  dgetri_(&n, NULL, &ldy, NULL, &size, &lwork, &info);
  if (info != 0)
    return info;

  // Allocate the workspace
  lwork = (int)size;
  double * work;
  if ((work = malloc((size_t)lwork * sizeof(double))) == NULL) {
    free(ipiv);
    return ENOMEM;
  }

  // Calculate the inverse
  dgetri_(&n, Y, &ldy, ipiv, work, &lwork, &info);

  free(ipiv);
  free(work);

  return info;
}

static mcmc_error calculate_equilibrium_state_occupancies(unsigned int n, const double * restrict Q, size_t ldq, double * restrict eqstates, ion_workspace * workspace) {
  double * S = workspace->X;
  size_t ldx = workspace->ldx;

  /* The original (Matlab) code does:
   *   u = ones(1, nstates);
   *   S = [ Q u' ];
   *   eqStates = u / (S * S');
   * This is equivalent to:
   *   S = ones(nstates);
   *   eq_states = ones(1, nstates);
   *   S += Q * Q';             // BLAS DGEMM
   *   eqstates *= inv(S);      // BLAS DGESV
   */
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++)
      S[j * lds + i] = 1.0;
    eqstates[j] = 1.0;
  }

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n,
              1.0, Q, ldq, Q, ldq, 1.0, S, lds);

  feclearexcept(FE_ALL_EXCEPT);
  cblas_dgesv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n,
              S, lds, eqstates, 1);

  // Check for numerical instability in the inverse
  if (fetestexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID) != 0)
    return MCMC_ERROR_NUMERICAL_INSTABILITY;

  return MCMC_SUCCESS;
}

static int calculate_specmat_eigenvectors(unsigned int n, const double * restrict Q, size_t ldq, double * restrict v, double * specMat, size_t lds) {
  double * X, * X_inv;
  if ((X = malloc(n * (ldx = ((n + 1u) & ~1u)) * sizeof(double))) == NULL)
    return ENOMEM;
  if (X_inv = malloc(n * ldx * sizeof(double))) == NULL) {
    free(X);
    return ENOMEM;
  }

  int error;
  if ((error = eigenvectors(n, Q, ldq, X, ldx, v)) != 0) {
    free(X);
    free(X_inv);
    return error;
  }
  if ((error = inverse(n, X, ldx, X_inv, ldx)) != 0) {
    free(X);
    free(X_inv);
    return error;
  }

  for (int j = 0; j < n; j++)
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, &X[j * ldx], ldx, &X_inv[j], ldx, &specMat[j * n * lds], lds);
  
  free(X);
  free(X_inv);
  
  return 0;
}
