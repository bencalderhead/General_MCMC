#ifndef UNIFORM_H
#define UNIFORM_H

#include <math.h>
#include "rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Uniform distribution over the open range (lower, upper).
 *
 * Parameters:
 *   lower < upper
 *
 * http://en.wikipedia.org/wiki/Uniform_distribution
 */
typedef struct {
  double lower, upper;
} uniform;

/**
 * Generates a uniformly distributed double precision floating point variable.
 *
 * @param r  a random number generator
 * @param u  distribution parameters
 *
 * @return a uniformly distributed double precision floating point variable.
 */
static inline double uniform_rand(const rng * restrict r, const uniform * restrict u) {
  return u->lower + rng_get_double(r) * (u->upper - u->lower);
}

/**
 * Calculates the value of the probability density function for the specified
 * uniform distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param u  distribution parameters
 *
 * @return the value of the PDF at point x.
 */
static inline double uniform_pdf(double x, const uniform * u) {
  if (x <= u->lower || x >= u->upper)
    return 0.0;
  return 1.0 / (u->upper - u->lower);
}

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified uniform distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param u  distribution parameters
 *
 * @return the value of the first derivative of the PDF at point x.
 */
static inline double uniform_1st_order_pdf(double x, const uniform * u) {
  if (x <= u->lower || x >= u->upper)
    return -HUGE_VAL;
  return 0.0;
}

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified uniform distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param u  distribution parameters
 * @return the value of the second derivative of the PDF at point x.
 */
static inline double uniform_2nd_order_pdf(double x, const uniform * u) {
  if (x <= u->lower || x >= u->upper)
    return -HUGE_VAL;
  return 0.0;
}

#ifdef __cplusplus
}
#endif

#endif
