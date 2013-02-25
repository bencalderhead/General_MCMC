#include "distributions/normal.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <errno.h>

#define RSQRT_2PI

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * Normal distribution.
 *
 * Parameters:
 *   mean
 *   variance > 0
 *
 * http://en.wikipedia.org/wiki/Normal_distribution
 */
struct __normal_distribution_st {
  double mean, variance;
};

/**
 * Creates a new normal distribution.
 *
 * @param n         the newly created normal distribution is returned through
 *                    this pointer
 * @param mean      the mean
 * @param variance  the variance
 *
 * @return 0 on success, EINVAL if variance <= 0, ENOMEM if there is not enough
 *  memory available to create another normal distribution.
 */
int normal_create(normal_distribution * n, double mean, double variance) {
  if (variance <= 0.0)
    return EINVAL;

  if ((*n = malloc(sizeof(struct __normal_distribution_st))) == NULL)
    return ENOMEM;

  (*n)->mean = mean;
  (*n)->variance = variance;

  return 0;
}

/**
 * Destroys a normal distribution, freeing any memory allocated.
 *
 * @param n  the normal distribution to destroy
 */
void normal_destroy(normal_distribution n) {
  free(n);
}

double normal_get_mean(const normal_distribution n) { return n->mean; }
double normal_get_variance(const normal_distribution n) { return n->variance; }

int normal_set_mean(normal_distribution n, double mean) {
  n->mean = mean;
  return 0;
}

int normal_set_variance(normal_distribution n, double variance) {
  if (variance <= 0.0)
    return EINVAL;
  n->variance = variance;
  return 0;
}

/**
 * Generates a normally distributed double precision floating point variable.
 * If the distribution is NULL then it is assumed to be N(0, 1).
 *
 * @param n      the normal distribution (may be NULL)
 * @param rand   a function generating double precision floating point variables
 *                 uniformly distributed over [0, 1).
 * @param state  state for the rng function
 *
 * @return a normally distributed double precision floating point variable.
 */
double normal_rand(const normal_distribution n, uint64_t (*rand)(const void *), const void * state) {
  // Box-Muller Transform
  static double next;
  static bool hasNext = false;

  if (hasNext) {
    hasNext = false;
    return (n == NULL) ? next : n->mean + next * n->variance;
  }

  double u0 = 1.0 - ((double)rand(state) / ((double)UINT64_MAX + 1.0));
  double u1 = 1.0 - ((double)rand(state) / ((double)UINT64_MAX + 1.0));
  double r = sqrt(-2.0 * log(u0));
  double theta = 2.0 * M_PI * u1;

  next = r * sin(theta);
  hasNext = true;

  double x = r * cos(theta);
  return (n == NULL) ? x : n->mean + x * n->variance;
}

/**
 * Calculates the log of the value of the probability density function for the
 * specified normal distribution at point x.  If the distribution is NULL then
 * it is assumed to be N(0, 1).
 *
 * @param n  the normal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the log of the value of the PDF at point x.
 */
double normal_log_pdf(const normal_distribution n, double x) {
  if (n == NULL)
    return -log(sqrt(2.0 * M_PI)) - ((x * x) / 2.0);
  return -log(n->variance * sqrt(2.0 * M_PI)) - (((x - n->mean) * (x - n->mean)) / (2.0 * n->variance * n->variance));
}

/**
 * Calculates the value of the probability density function for the specified
 * normal distribution at point x.  If the distribution is NULL then it is
 * assumed to be N(0, 1).
 *
 * @param n  the normal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the PDF at point x.
 */
double normal_pdf(const normal_distribution n, double x) {
  return (1.0 / (n->variance * 2.0 * M_PI)) * exp(-((x - n->mean) * (x - n->mean)) / (2.0 * n->variance * n->variance));
}

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified normal distribution at point x.
 *
 * @param n  the normal distribution
 * @param x  the point to evaluate the PDF at
 * @return the value of the first derivative of the PDF at point x.
 */
double normal_1st_order_pdf(const normal_distribution n, double x) {
  return -(x - n->mean) / (n->variance * n->variance);
}

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified normal distribution at point x.
 *
 * @param n  the normal distribution
 * @param x  the point to evaluate the PDF at
 * @return the value of the second derivative of the PDF at point x.
 */
double normal_2nd_order_pdf(const normal_distribution n, double x) {
  (void)x;
  return -1.0 / (n->variance * n->variance);
}
