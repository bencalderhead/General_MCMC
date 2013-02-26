#ifndef NORMAL_H
#define NORMAL_H

#include <math.h>
#include "rng.h"

// M_PI is not defined in C99
// M_PI always appears as sqrt(2.0 * M_PI) so define that here instead
// This is also used by the lognormal and gamma distributions
#define M_SQRT2PI 2.50662827463100050241

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Normal distribution.
 *
 * Parameters:
 *   mean
 *   standard deviation > 0
 *
 * http://en.wikipedia.org/wiki/Normal_distribution
 */
typedef struct {
  double mean, stddev;
} normal;

/**
 * Generates a normally distributed double precision floating point variable.
 *
 * @param r  a random number generator
 * @param n  distribution parameters
 *
 * @return a normally distributed double precision floating point variable.
 */
static inline double normal_rand(const rng * r, const normal * n) {
  // Numerical Recipes in C++, 3rd edition (section 7.3, page 369)

  double u, v, q;
  do {
    u = rng_get_double(r);
    v = 1.7156 * (rng_get_double(r) - 0.5);
    double x = u - 0.449871;
    double y = fabs(v) + 0.386595;
    q = x * x + y * (0.19600 * y - 0.25472 * x);
  } while (q > 0.27597 && (q > 0.27846 || v * v > -4.0 * log(u) * u * u));

  return n->mean + n->stddev * v / u;
}

/**
 * Calculates the log of the value of the probability density function for the
 * specified normal distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param n  distribution parameters
 *
 * @return the log of the value of the PDF at point x.
 */
static inline double normal_log_pdf(double x, const normal * n) {
  return -log(n->stddev) - log(M_SQRT2PI) -
  ((x - n->mean) * (x - n->mean)) / (2.0 * n->stddev * n->stddev);
}

/**
 * Calculates the value of the probability density function for the specified
 * normal distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param n  distribution parameters
 *
 * @return the value of the PDF at point x.
 */
static inline double normal_pdf(double x, const normal * n) {
  return 1.0 / (n->stddev * M_SQRT2PI) *
  exp(-((x - n->mean) * (x - n->mean)) / (2.0 * n->stddev * n->stddev));
}

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified normal distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param n  distribution parameters
 *
 * @return the value of the first derivative of the PDF at point x.
 */
static inline double normal_1st_order_pdf(double x, const normal * n) {
  return -(x - n->mean) / (n->stddev * n->stddev);
}

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified normal distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param n  distribution parameters
 *
 * @return the value of the second derivative of the PDF at point x.
 */
static inline double normal_2nd_order_pdf(double x, const normal * n) {
  (void)x;
  return -1.0 / (n->stddev * n->stddev);
}

#ifdef __cplusplus
}
#endif

#endif
