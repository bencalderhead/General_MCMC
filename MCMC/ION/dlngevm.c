#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

static size_t max(size_t a, size_t b) { return (a > b) ? a : b; }

#ifdef MATLAB_MEX_FILE
#include "mex.h"

static void xerbla(const char * func, int info) {
  char buf[256];
  snprintf(buf, 256, "On entry to %s parameter %d had an invalid value", func, info);
  mexErrMsgTxt(buf);
}

#define calloc(n, size) mxCalloc(n, size)
#define free(ptr) mxFree(ptr)

#else

static void xerbla(const char * func, int info) {
  fprintf(stderr, "On entry to %s parameter %d had an invalid value", func, info);
}

#endif

#define XERBLA(info) xerbla(__func__, info)

/**
 * Compares x to y as doubles.  Used to sort double values into descending
 * order.
 *
 * @param x,y  two doubles to compare
 * @return -1 if x > y, 1 if x < y, 0 otherwise.
 */
static int comparator(const void * x, const void * y) {
  double a, b;

  a = *(double*)x;
  b = *(double*)y;

  return (a < b) ? 1 : ((a == b) ? 0 : -1);
}

/**
 * Performs:
 *      y = x * A
 * in log space.
 *
 * @param m     number of rows in A/elements in x
 * @param n     number of columns in A/elements in y
 * @param x     array of elements in x in log space
 * @param incx  spacing between elements of x
 * @param A     array of elements in A in column-major layout
 * @param lda   spacing between columns of A
 * @param y     array of elements in y (uninitialised)
 * @param incy  spacing between elements of y
 */
void dlngevm(size_t m, size_t n, const double * x, size_t incx, const double * A, size_t lda, double * y, size_t incy) {
  size_t i, j;
  int info;
  double * temp, sum;

  info = 0;
  if (incx < 1)
    info = 4;
  else if (lda < max(1, m))
    info = 6;
  else if (incy < 1)
    info = 8;
  if (info != 0) {
    XERBLA(info);
    return;
  }

  if (m == 0 || n == 0) return;

  temp = calloc(m, sizeof(double));     /* Reset to zeros */
  if (temp == NULL) {
    fputs("Error allocating temporary array\n", stderr);
    return;
  }

  for (j = 0; j < n; j++) {
    /* Multiply row by column in log - i.e. sum! */
    for (i = 0; i < m; i++)
      temp[i] = x[i * incx] + A[j * lda + i];

    if (m > 1) {
      /* Now add them all together in log to get jth element */
      qsort(temp, m, sizeof(double), comparator); /* Sort with biggest first */

      /* Check the first values is not -inf */
      if (temp[0] != -HUGE_VAL) {
        /* Then do normal exponential sum */
        sum = 0.0;
        for (i = 1; i < m; i++)
          sum += exp(temp[i] - temp[0]);
        y[j * incy] = temp[0] + log(1.0 + sum);
      }
      else
        y[j * incy] = -HUGE_VAL;
    }
    else
      /* Only one value so no need to sum! */
      y[j * incy] = temp[0];
  }
  
  free(temp);
}

#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) {
  const double * x, * A;
  double * y;
  size_t m, n;

  if (nlhs != 1)
    mexErrMsgTxt("One output argument required");
  if (nrhs != 2)
    mexErrMsgTxt("Two input arguments required");

  x = mxGetPr(prhs[0]);
  A = mxGetPr(prhs[1]);
  m = mxGetN(prhs[0]);
  n = mxGetN(prhs[1]);
  y = mxGetPr(plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL));

  if (m != mxGetM(prhs[1]))
    mexErrMsgTxt("Inner dimensions do not match");

  dlngevm(m, n, x, 1, A, m, y, 1);
}

#else
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>

static void dgevm(size_t, size_t, const double *, size_t, const double *, size_t, double *, size_t);

int main(int argc, char * argv[]) {
  double * x, * lnx, * y, * lny, * A, * lnA, diff, d, time;
  size_t m, n, incx, incy, lda, i, j;
  struct timeval start, stop;

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <m> <n>\nwhere:\n  m and n         are the sizes of the matrices\n", argv[0]);
    return 1;
  }

  if (sscanf(argv[1], "%lu", &m) != 1) {
    fprintf(stderr, "Unable to parse number from '%s'\n", argv[1]);
    return 1;
  }

  if (sscanf(argv[2], "%lu", &n) != 1) {
    fprintf(stderr, "Unable to parse number from '%s'\n", argv[2]);
    return 2;
  }

  srand(0);

  if ((x = malloc(m * (incx = 1) * sizeof(double))) == NULL) {
    fputs("Unable to allocate x\n", stderr);
    return -1;
  }

  if ((A = malloc((lda = (m + 3u) & ~3u) * n * sizeof(double))) == NULL) {
    fputs("Unable to allocate A\n", stderr);
    return -2;
  }

  if ((y = malloc(n * (incy = 1) * sizeof(double))) == NULL) {
    fputs("Unable to allocate y\n", stderr);
    return -3;
  }

  if ((lnx = malloc(m * (incx = 1) * sizeof(double))) == NULL) {
    fputs("Unable to allocate lnx\n", stderr);
    return -4;
  }

  if ((lnA = malloc((lda = (m + 3u) & ~3u) * n * sizeof(double))) == NULL) {
    fputs("Unable to allocate lnA\n", stderr);
    return -5;
  }

  if ((lny = malloc(n * (incy = 1) * sizeof(double))) == NULL) {
    fputs("Unable to allocate lny\n", stderr);
    return -6;
  }

  for (i = 0; i < m; i++)
    lnx[i * incx] = log(x[i * incx] = (double)rand() / (double)RAND_MAX);

  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      lnA[j * lda + i] = log(A[j * lda + i] = (double)rand() / (double)RAND_MAX);
  }

  dgevm(m, n, x, incx, A, lda, y, incy);
  dlngevm(m, n, lnx, incx, lnA, lda, lny, incy);

  fputs("x = \n", stdout);
  for (i = 0; i < m; i++)
    fprintf(stdout, "%15.6f", x[i * incx]);
  fputs("\n", stdout);

  fputs("A = \n", stdout);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      fprintf(stdout, "%15.6f", A[j * lda + i]);
    fputs("\n", stdout);
  }

  fputs("y = \n", stdout);
  for (j = 0; j < n; j++)
    fprintf(stdout, "%15.6f\n", y[j * incy]);

  fputs("log(y) = \n", stdout);
  for (j = 0; j < n; j++)
    fprintf(stdout, "%15.6f\n", log(y[j * incy]));

  fputs("lnx = \n", stdout);
  for (i = 0; i < m; i++)
    fprintf(stdout, "%15.6f", lnx[i * incx]);
  fputs("\n", stdout);

  fputs("lnA = \n", stdout);
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
      fprintf(stdout, "%15.6f", lnA[j * lda + i]);
    fputs("\n", stdout);
  }

  fputs("lny = \n", stdout);
  for (j = 0; j < n; j++)
    fprintf(stdout, "%15.6f\n", lny[j * incy]);

  diff = 0.0;
  for (j = 0; j < n; j++) {
    d = fabs(log(y[j * incy]) - lny[j * incy]);
    if (d > diff)
      diff = d;
  }
  fprintf(stdout, "Error = %.3e\n", diff);
  if (diff <= (3 * m - 1) * DBL_EPSILON)
    fputs("PASSED!\n", stderr);
  else
    fputs("FAILED!\n", stderr);

  if (gettimeofday(&start, NULL) != 0) {
    fputs("Get start time failed\n", stderr);
    return -7;
  }
  for (i = 0; i < 20; i++)
    dlngevm(m, n, lnx, incx, lnA, lda, lny, incy);
  if (gettimeofday(&stop, NULL) != 0) {
    fputs("Get stop time failed\n", stderr);
    return -8;
  }

  time = ((double)(stop.tv_sec - start.tv_sec) + (double)(stop.tv_usec - start.tv_usec) * 1.e-6) / 20.0;
  fprintf(stdout, "Time = %.3es\n", time);
  fprintf(stdout, "%.3eGFlops/s\n", (n * (3 * m - 1) * 1.e-9) / time);

  free(x);
  free(y);
  free(A);
  free(lnx);
  free(lny);
  free(lnA);

  return 0;
}

static void dgevm(size_t m, size_t n, const double * x, size_t incx, const double * A, size_t lda, double * y, size_t incy) {
  size_t i, j;
  int info;
  double temp;

  info = 0;
  if (incx < 1)
    info = 4;
  else if (lda < max(1, m))
    info = 6;
  else if (incy < 1)
    info = 8;
  if (info != 0) {
    XERBLA(info);
    return;
  }

  if (m == 0 || n == 0) return;

  for (j = 0; j < n; j++) {
    temp = 0.0;
    for (i = 0; i < m; i++)
      temp += x[i * incx] * A[j * lda + i]; /* So much simpler in linear space! */
    y[j * incy] = temp;
  }

}

#endif
