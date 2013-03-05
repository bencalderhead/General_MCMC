#include <stdlib.h>
#include <complex.h>
#include <errno.h>

extern void dgeev_(const char *, const char *, const int *,
                   double *, const int *,
                   double *, double *,
                   double *, const int *, double *, const int *, double *,
                   const int *, int *);

/**
 * Calculates the eigenvectors and eigenvalues of an n-by-n real nonsymmetric
 * matrix A using the QR decomposition.
 *
 * The computed eigenvectors are normalized to have Euclidean norm equal to 1
 * and largest component real.
 *
 * @param n    the size of the matrix A.
 * @param A    a pointer to a double precision array of size n * ldx containing
 *               the elements of A in column major layout. A is overwritten on
 *               output.
 * @param lda  the leading dimension of the array A.
 * @param B    a pointer to a double precision complex array of size n * ldb.
 *               The eigenvectors of A are stored column by column in B.
 * @param ldb  the leading dimension of the array B.
 * @param x    a pointer to a double precision complex array of size n.  The
 *               eigenvalues of A are stored in x.
 *
 * @return 0 on success,
 *         EINVAL if the parameters were invalid,
 *         ENOMEM if the workspace could not be allocated,
 *         ERANGE if the eigenvalues or eigenvectors could not be calculated.
 */
int eigen(int n, double * restrict A, int lda,
          double complex * restrict B, int ldb, double complex * x) {
  int info = 0;

  /* Parameter checking */
  if (n < 0)
    info = -1;
  else if (lda < n)
    info = -3;
  else if (ldb < n)
    info = -5;

  if (info != 0)
    return info;

  /* Quick return if possible */
  if (n == 0)
    return info;

  // Calling DGEEV with lwork = -1 causes it to return the optimal workspace size in work[0]
  double size;
  int lwork = -1;
  dgeev_("N", "V", &n, NULL, &lda, NULL, NULL, NULL, &ldb, NULL, &ldb, &size, &lwork, &info);

  double * wr, * wi, * work;
  if ((wr = malloc((size_t)n * sizeof(double))) == NULL)
    return ENOMEM;
  if ((wi = 
}


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
