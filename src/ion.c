#include "ion.h"
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fenv.h>
#include <errno.h

#define ERROR_CHECK(call) \
  do { \
    int error = (call); \
    if (error != 0) \
      return error; \
  } while (false)

static inline int vector_alloc(int n, void ** x, int incx, size_t size) {
  *x = malloc((size_t)(n * incx) * size);
  return (*x == NULL) ? ENOMEM : 0;
}

#define vector_init(n, alpha, x, incx) \
  do { \
    for (int i = 0; i < (n); i++) \
      (x)[i * (incx)] = (alpha); \
  } while (false)

#define vector_copy(n, i, x, incx, y, incy) \
  do { \
    if ((i) == NULL) { \
      for (int ii = 0; ii < (n); ii++) \
        (y)[ii * (incy)] = (x)[ii * (incx)]; \
    } \
    else { \
      for (int ii = 0; ii < (n); ii++) \
        (y)[ii * (incy)] = (x)[(i)[ii] * (incx)]; \
    } \
  } while (false)

#define vector_log(n, x, incx, y, incy) \
  do { \
    for (int ii = 0; ii < (n); ii++) \
      (y)[ii * (incy)] = log((x)[ii * (incx)]); \
  } while (false)

static inline int matrix_alloc(int m, int n, void ** A, int * lda, size_t size) {
  const unsigned int align = 16u / size;
  *lda = (int)(((unsigned int)m + align - 1u) & ~(align - 1u));
  *A = malloc((size_t)(n * *lda) * size);
  return (*A == NULL) ? OUT_OF_MEMORY : 0;
}

#define matrix_init(m, n, alpha, A, lda) \
  do { \
    for (int j = 0; j < (n); j++) { \
      for (int i = 0; i < (m); i++) \
        (A)[j * (lda) + i] = (alpha); \
    } \
  } while (false)

#define matrix_copy(m, n, i, j, A, lda, B, ldb) \
  do { \
    if ((i) == NULL || (j) == NULL) { \
      for (int jj = 0; jj < (n); jj++) { \
        for (int ii = 0; ii < (n); ii++) \
          (B)[jj * (ldb) + ii] = (A)[jj * (lda) + ii]; \
      } \
    } \
    else { \
      for (int jj = 0; jj < (n); jj++) { \
        for (int ii = 0; ii < (n); ii++) \
          (B)[jj * (ldb) + ii] = (A)[(j)[jj] * (lda) + (i)[ii]]; \
      } \
    } \
  } while (false)

#define diag(n, A, lda, x, incx) \
  do { \
    for (int i = 0; i < (n); i++) \
      (x)[i * (incx)] = (A)[i * (lda) + i]; \
  } while (false)

static int eig(int, const double *, int, double complex *, int, double complex *, int);
static int inv(int, const double *, int, double *, int);

static void outer_product(int, const double *, int, const double *, int, double *, int);

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
  double * ll;          // log likelihood
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
  int ldq;
  ERROR_CHECK(matrix_alloc(model->nStates, model->nStates, (void **)&Q, &ldq, sizeof(double)));
  model->calculate_Q_matrix(chain->params, Q, ldq);

  // Split up the Q matrix into its component matrices  
  double * Q_FF, * Q_FA, * Q_AF, * Q_AA;
  int ldqff, ldqfa, ldqaf. ldqaa;
  ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nClosedStates, (void **)&Q_FF, &ldqff, sizeof(double)));
  ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nOpenStates,   (void **)&Q_FA, &ldqfa, sizeof(double)));
  ERROR_CHECK(matrix_alloc(model->nOpenStates,   model->nClosedStates, (void **)&Q_AF, &ldqaf, sizeof(double)));
  ERROR_CHECK(matrix_alloc(model->nOpenStates,   model->nOpenStates,   (void **)&Q_AA, &ldqaa, sizeof(double)));

  matrix_copy(model->nClosedStates, model->nClosedStates, Q, ldq, Q_FF, ldqff);
  matrix_copy(model->nClosedStates, model->nOpenStates,   Q, ldq, Q_FA, ldqfa);
  matrix_copy(model->nOpenStates,   model->nClosedStates, Q, ldq, Q_AF, ldqaf);
  matrix_copy(model->nOpenStates,   model->nOpenStates, Q, ldq, Q_AA, ldqaa);

  // Calculate equilibrium state occupancies
  // Original code does:
  //  u = ones(1, numstates);
  //  S = [ Q u' ];
  //  EqStates = u / (S * S');
  // This is equivalent to:
  //  u = ones(1, numstates);
  //  eqstates = u / (Q * Q' + ones(numstates));

  double * S, * eqStates;
  int lds;
  ERROR_CHECK(matrix_alloc(model->nStates, model->nStates, (void **)&S, &lds, sizeof(double)));
  ERROR_CHECK(vector_alloc(model->nStates, (void **)&eqStates, 1, sizeof(double)));
  
  // Initialise S to all ones (S = ones(numstates))
  matrix_init(model->nStates, model->nStates, 1.0, S, lds);
  // Initialise eqStates to all ones (eqStates = ones(1, numstates))
  vector_init(model->nStates, 1.0, eqStates, 1);

  // S += Q * Q' (S = Q * Q' + ones(numstates))
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, model->nStates, model->nStates, model->nStates, 1.0, Q, ldq, Q, ldq, 1.0, S, lds);
  free(Q);
  // eqStates /= S (eqStates = ones(1, numstates) / (Q * Q' + ones(numstates)))
  feclearexcept(FE_DIVBYZERO);
  cblas_dgesv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, model->nStates, S, lds, eqStates, 1);      // Numerical instability may occur in here
  free(S);
  if (fegetexcept(FE_DIVBYZERO) != 0)
    return ION_CHANNEL_NUMERICAL_INSTABILITY;

  double * eqStates_F, * eqStates_A;
  ERROR_CHECK(vector_alloc(model->nClosedStates, &eqStates_F, 1, sizeof(double)));
  ERROR_CHECK(vector_alloc(model->nOpenStates,   &eqStates_A, 1, sizeof(double)));
  vector_log(model->nClosedStates, model->closedStates, eqStates, 1, eqStates_F, 1);
  vector_log(model->nOpenStates, model->openStates, eqStates, 1, eqStates_A, 1);

  // Calculate spectral matrices and eigenvectors of current Q_FF
  double * X, * V, * Y, * V_Q_FF, ** SpecMat_Q_FF, * V_Q_AA, ** SpecMat_Q_AA;
  int ldx, ldv, ldy, * ldspecMat_Q_FF, * ldspecMat_Q_AA;
  ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nClosedStates, (void **)&X, &ldx, sizeof(double)));
  ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nClosedStates, (void **)&V, &ldv, sizeof(double)));
  ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nClosedStates, (void **)&Y, &ldy, sizeof(double)));
  ERROR_CHECK(vector_alloc(model->nClosedStates, (void **)&V_Q_FF, 1, sizeof(double)));
  ERROR_CHECK(vector_alloc(model->nClosedStates, (void **)&SpecMat_Q_FF, 1, sizeof(double *)));
  ERROR_CHECK(vector_alloc(model->nClosedStates, (void **)&ldspecMat_Q_FF, 1, sizeof(int)));

  ERROR_CHECK(eig(model->nClosedStates, Q_FF, ldqff, X, ldx, V, ldv));
  diag(model->nClosedStates, V, ldv, V_Q_FF, 1);        // Eigenvectors
  ERROR_CHECK(inv(model->nClosedStates, X, ldx, Y, ldy));
  for (int j = 0; j < model->nClosedStates; j++) {
    ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nClosedStates, (void **)&SpecMat_Q_FF[j], &ldspecMat_Q_FF[j], sizeof(double)));
    outer_product(model->nClosedStates, &X[j * ldx], 1, &Y[j], ldy, &SpecMat_Q_FF[j], ldspecMat_Q_FF[j]);    // Calculate spectral matrices
  }

  free(X);
  free(V);
  free(Y);

  // Calculate spectral matrices and eigenvectors of current Q_AA
  ERROR_CHECK(matrix_alloc(model->nOpenStates, model->nOpenStates, (void **)&X, &ldx, sizeof(double)));
  ERROR_CHECK(matrix_alloc(model->nOpenStates, model->nOpenStates, (void **)&V, &ldv, sizeof(double)));
  ERROR_CHECK(matrix_alloc(model->nOpenStates, model->nOpenStates, (void **)&Y, &ldy, sizeof(double)));
  ERROR_CHECK(vector_alloc(model->nOpenStates, (void **)&V_Q_AA, 1, sizeof(double)));
  ERROR_CHECK(vector_alloc(model->nOpenStates, (void **)&SpecMat_Q_AA, 1, sizeof(double *)));
  ERROR_CHECK(vector_alloc(model->nOpenStates, (void **)&ldspecMat_Q_AA, 1, sizeof(int)));

  ERROR_CHECK(eig(model->nOpenStates, Q_AA, ldqaa, X, ldx, V, ldv));
  diag(model->nOpenStates, V, ldv, V_Q_FF, 1);        // Eigenvectors
  ERROR_CHECK(inv(model->nOpenStates, X, ldx, Y, ldy));
  for (int j = 0; j < model->nOpenStates; j++) {
    ERROR_CHECK(matrix_alloc(model->nOpenStates, model->nOpenStates, (void **)&SpecMat_Q_AA[j], &ldspecMat_Q_AA[j], sizeof(double)));
    outer_product(model->nOpenStates, &X[j * ldx], 1, &Y[j], ldy, &SpecMat_Q_AA[j], ldspecMat_Q_AA[j]);    // Calculate spectral matrices
  }

  free(X);
  free(V);
  free(Y);


  // Calculate initial vectors for current state
  double * L =  (model->data[0] == 0) ? eqStates_F : eqStates_A;


  // Do in a slow loop to begin with
  for (int i = 0; i < model->nTimePoints - 1; i++) {
    if (mode->data[i] == 0) {   // Currently in closed state (denoted F in literature)

      // Moving to open state
      double sojourn = model->timePoints[i + 1] - model->timePoints[i]; // Get time interval to next move

      // Calculate idealised transition probability from closed to open

      double * G_FA;
      int ldgfa;
      ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nClosedStates, (void **)&G_FA, &ldgfa, sizeof(double)));
      matrix_init(model->nClosedStates, model->nClosedStates, 0.0, G_FA, ldgfa);

      for (int j = 0; j < model->nClosedStates; j++) {
        double alpha = exp(V_Q_FF[j] * sojourn);
        for (int k = 0; k < model->nClosedStates; k++)
          cblas_daxpy(model->nClosedStates, alpha, &SpecMat_Q_FF[j][k * ldspecMat_Q_FF[j]], 1, &G_FA[k * ldgfa], 1);
      }

      double * T;
      int ldt;
      ERROR_CHECK(matrix_alloc(model->nClosedStates, model->nClosedStates, (void **)&T, &ldt, sizeof(double)));
      matrix_copy(model-nClosedStates, model-nClosedStates, G_FA, ldgfa, T, ldt);
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, model-nClosedStates, model-nClosedStates, model-nClosedStates, 1.0, T, ldt, Q_FA, ldqfa, 1.0, G_FA, ldqfa);
      free(T);

      // Logarithmic update
      // Calculate log(L*G_FA) in terms of LL = log(L) and log(G_FA)
      Log_G_FA = log(G_FA);

      // Do element-wise logarithmic calculation
      // log (vector * matrix) in terms of log(vector) and log(matrix)
      LL = dlngevm(LL, Log_G_FA);
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

  lwork = (int)size;
  double * work;
  if ((work = malloc((size_t)lwork * sizeof(double))) == NULL)
    return ENOMEM;

  // Calculate the inverse
  dgetri_(&n, Y, &ldy, ipiv, work, &lwork, &info);

  free(ipiv);
  free(work);

  return info;
}

static void outer_product(int n, const double * x, int incx, const double * y, int incy, double * A, int lda) {
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++)
      A[j * lda + i] = x[i * incx] * y[j * incy];
  }
}
