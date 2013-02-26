#ifndef LOGNORMAL_H
#define LOGNORMAL_H

#include <math.h>
#include "rng.h"
#include "normal.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * LogNormal distribution.
 *
 * Parameters:
 *   shape > 0
 *   log-scale
 *
 * http://en.wikipedia.org/wiki/Lognormal_distribution
 */
typedef struct {
  double mean, stddev;
} lognormal;

/**
 * Generates a lognormal distributed double precision floating point variable.
 *
 * @param r  a random number generator
 * @param l  distribution parameters
 *
 * @return a lognormal distributed double precision floating point variable.
 */
static inline double lognormal_rand(const rng * r, const lognormal * l) {
  const normal n = { 0.0, 1.0 };
  return exp(l->mean + l->stddev * normal_rand(r, &n));
}

/**
 * Calculates the value of the probability density function for the specified
 * lognormal distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param l  distribution parameters
 *
 * @return the value of the PDF at point x.
 */
static inline double lognormal_pdf(double x, const lognormal * l) {
  return (1.0 / (x * l->stddev * M_SQRT2PI)) *
  exp(-((log(x) - l->mean) * (log(x) - l->mean)) / (2.0 * l->stddev * l->stddev));
}

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified lognormal distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param l  distribution parameters
 *
 * @return the value of the first derivative of the PDF at point x.
 */
static inline double lognormal_1st_order_pdf(double x, const lognormal * l) {
  return -(1.0 / x) - ((log(x) - l->mean) / (l->stddev * l->stddev)) * (1.0 / x);
}

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified lognormal distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param l  distribution parameters
 *
 * @return the value of the second derivative of the PDF at point x.
 */
static inline double lognormal_2nd_order_pdf(double x, const lognormal * l) {
  return (1.0 / (x * x)) + ((log(x) - l->mean) / (l->stddev * l->stddev)) * (1.0 / (x * x)) - (1.0 / (l->stddev * l->stddev)) * (1.0 / (x * x));
}

#ifdef __cplusplus
}
#endif

#endif
