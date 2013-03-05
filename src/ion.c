#include "ion.h"
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fenv.h>
#include <errno.h>

struct __ion_workspace_st {
  double * Q, * eq_states, * S;
  size_t ldq;

};

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
static void eigen(int, double *, int, double complex *, int, double complex *);
static int inv(int, const double *, int, double *, int);

int ion_model_create(ion_model * restrict model, unsigned int open, unsigned int closed, void (*update_q)(const double *, double *, size_t)) {
  int error = mcmc_model_create(&model->m);
  if (error != 0)
    return error;

  model->open = open;
  model->closed = closed;
  model->update_q = update_q;

  // Work out the workspace size
  size_t w = 0;

  // n is the total number of states
  unsigned int n = open + closed;
  size_t ld = ((n + 1u) & ~1u) * sizeof(double);        // n rounded up to the alignment

  w += 2 * n * ld + ld;

  ld = ((closed + 1u) & ~1u) * sizeof(double);
  w +=

  if ((model->workspace = malloc(sizeof(ion_workspace))) == NULL)
    return ENOMEM;

  return 0;
}

void ion_model_destroy(ion_model * model) {
  mcmc_model_destroy(&model->m);
  return 0;
}

int ion_evaluate_mh(const ion_model * restrict model, markov_chain * restrict chain) {
  int error;

  // Calculate the prior for each parameter
  for (int i = 0; i < model->n; i++) {
    double x = prior_evaluate(&model->priors[i], chain->params[i]);
    if (isnan(x) || islessequal(x, 0.0))
      return EDOM;
    chain->logprior[i] = log(x);
  }


  // Calculate the likelihood of the ion channel data

  unsigned int states = model->open + model->closed;

  // Allocate memory for the Q matrix and equilibrium states
  double * Q = model->workspace;                // states * ldq
  double * eq_states = &Q[states * ldq];        // states
  size_t ldq = (states + 1u) & ~1u;

  // Calculate the Q matrix
  model->update_q(chain->params, Q, ldq);

  // Calculate equilibrium state occupancies
  double * S = &eq_states[ldq];                 // states * lds
  size_t lds = ldq;

  for (int j = 0; j < states; j++) {
    for (int i = 0; i < states; i++)
      S[j * lds + i] = 1.0;
    eq_states[j] = 1.0;
  }

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, states, states, states, 1.0, Q, ldq, Q, ldq, 1.0, S, lds);

  feclearexcept(FE_ALL_EXCEPT);
  cblas_dgesv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, states, S, lds, eq_states, 1);
  if (fetestexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID) != 0)
    return ERANGE;

  // Split up the Q matrix into its component matrices
  double * Q_FF = Q;
  double * Q_FA = &Q[model->open * ldq];
  double * Q_AF = &Q[model->open];
  double * Q_AA = &Q[model->open * ldq + model->open];

  // Split up the equilibrium states
  double * eq_states_F = eq_states;
  double * eq_states_A = &eq_states[model->open];


  // Calculate spectral matrices and eigenvectors of current Q_FF
  double * v_Q_FF = &eq_states[];
  double * specMat_Q_FF = malloc(model->closed * model->closed * ldsqff * sizeof(double));
  size_t ldsqff = (model->closed + 1u) & ~1u;

  size_t ldx = (model->closed + 1u) & ~1u;
  double * X = malloc(model->closed * ldx * sizeof(double));
  eig(model->closed, Q_FF, ldq, X, ldx, v_Q_FF);

  size_t ldy = (model->closed + 1u) & ~1u;
  double * Y = malloc(model->closed * ldy * sizeof(double));
  inv(model->closed, X, ldx, Y, ldy);

  for (int j = 0; j < n; j++)
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                model->closed, model->closed, model->closed,
                1.0, &X[j * ldx], ldx, Y[j], ldy,
                0.0, &specMat_Q_FF[j * model->closed * lds], ldsqff);

  free(X);
  free(Y);

  // Calculate spectral matrices and eigenvectors of current Q_AA
  size_t ldsqaa = (model->open + 1u) & ~1u;
  double complex * v_Q_AA = malloc(model->open * sizeof(double complex));
  double complex * specMat_Q_AA = malloc(model->open * model->open * ldsqaa * sizeof(double complex));

  ldx = (model->open + 1u) & ~1u;
  X = malloc(model->open * ldx * sizeof(double complex));
  eig(model->open, Q_AA, ldq, X, ldx, v_Q_AA);

  ldy = (model->open + 1u) & ~1u;
  Y = malloc(model->open * ldy * sizeof(double complex));
  inv(model->open, X, ldx, Y, ldy);

  for (int j = 0; j < n; j++)
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                model->open, model->open, model->open,
                1.0, &X[j * ldx], ldx, Y[j], ldy,
                0.0, &specMat_Q_AA[j * model->open * lds], ldsqaa);

  free(X);
  free(Y);


  // Calculate initial vectors for current state
  double * L =  (model->data[0] == 0) ? eqStates_F : eqStates_A;


  mcmc_model * m = model->m;
  // Do in a slow loop to begin with
  for (int i = 0; i < m->n_time_points - 1; i++) {
    double sojourn = m->time_points[i + 1] - m->time_points[i]; // Get time interval to next move

    if (m->data[i] == 0) {   // Currently in closed state (denoted F in literature)

      // Calculate idealised transition probability from closed to open

      size_t ldgfa = (model->closed + 1u) & ~1u;
      double complex * G_FA = malloc(model->closed * ldgfa * sizeof(double complex));

      for (int j = 0; j < model->closed; j++) {
        for (int i = 0; i < model->closed; i++)
          G_FA[j * ldgfa + i] = 0.0;
      }

      for (int k = 0; k < model->closed; k++) {
        double alpha = cexp(V_Q_FF[j] * sojourn);
        double complex * specMat = &specMat_Q_FF[k * model->closed * ldsqff];
        for (int j = 0; j < model->closed; j++) {
          for (int i = 0; i < model->closed; i++)
            G_FA[j * lda + i] += alpha * specMat[j * ldsqff + i];
        }
      }

      size_t ldt = (model->closed + 1u) & ~1u;
      double complex * T = malloc(model->closed * ldt * sizeof(double complex));
      for (int j = 0; j < model->closed; j++) {
        for (int i = 0; i < model->closed; i++)
          T[j * ldt + i] = G_FA[j * ldgfa + i];
      }
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, model->closed, model->closed, model->closed, 1.0, T, ldt, Q_FA, ldqfa, 1.0, G_FA, ldqfa);
      free(T);

      // Logarithmic update
      // Calculate log(L*G_FA) in terms of LL = log(L) and log(G_FA)
      // Do element-wise logarithmic calculation
      // log (vector * matrix) in terms of log(vector) and log(matrix)
      cblas_dlngevm(CblasTrans, 1.0, G_FA, ldgfa, LL, 1, 0.0, LL, 1);
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

static int calculate_equilibrium_state_occupancies(unsigned int n, const double * restrict Q, size_t ldq, double * restrict eq_states) {
  size_t lds = (n + 1u) & ~1u;
  double * S = malloc(n * lds * sizeof(double));
  if (S == NULL)
    return ENOMEM;

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++)
      S[j * lds + i] = 1.0;
    eq_states[j] = 1.0;
  }

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, Q, ldq, Q, ldq, 1.0, S, lds);

  feclearexcept(FE_ALL_EXCEPT);
  cblas_dgesv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, S, lds, eq_states, 1);

  free(S);

  if (fetestexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID) != 0)
    return ERANGE;

  return 0;
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
