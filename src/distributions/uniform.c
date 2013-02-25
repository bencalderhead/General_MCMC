#include "distributions/uniform.h"
#include <stdlib.h>
#include <math.h>
#include <errno.h>

/**
 * Uniform distribution.
 *
 * Parameters:
 *   lower < upper
 *
 * http://en.wikipedia.org/wiki/Uniform_distribution
 */
struct __uniform_distribution_st {
  double lower, upper;
};

/**
 * Creates a new uniform distribution over [lower, upper).
 *
 * @param u      the newly created uniform distribution is returned through this
 *                 pointer
 * @param lower  the lower bound
 * @param upper  the upper bound
 *
 * @return 0 on success, EINVAL if lower >= upper, ENOMEM if there is not enough
 *  memory available to create another uniform distribution.
 */
int uniform_create(uniform_distribution * u, double lower, double upper) {
  if (lower >= upper)
    return EINVAL;

  if ((*u = malloc(sizeof(struct __uniform_distribution_st))) == NULL)
    return ENOMEM;

  (*u)->lower = lower;
  (*u)->upper = upper;

  return 0;
}

/**
 * Destroys a uniform distribution, freeing any memory allocated.
 *
 * @param u  the uniform distribution to destroy
 */
void uniform_destroy(uniform_distribution u) {
  free(u);
}

double uniform_get_lower(const uniform_distribution u) { return u->lower; }
double uniform_get_upper(const uniform_distribution u) { return u->upper; }

int uniform_set_lower(uniform_distribution u, double lower) {
  if (lower >= u->upper)
    return EINVAL;
  u->lower = lower;
  return 0;
}

int uniform_set_upper(uniform_distribution u, double upper) {
  if (u->lower >= upper)
    return EINVAL;
  u->upper = upper;
  return 0;
}

/**
 * Generates a uniformly distributed double precision floating point variable.
 * If the distribution is NULL then it is assumed to be U[0, 1).
 *
 * @param u      the uniform distribution (may be NULL)
 * @param rand   a function generating double precision floating point variables
 *                 uniformly distributed over [0, 1)
 * @param state  state for the rng function
 *
 * @return a uniformly distributed double precision floating point variable.
 */
double uniform_rand(const uniform_distribution u, double (*rand)(const void *), const void * state) {
  double x = rand(state);
  return (u == NULL) ? x : u->lower + x * (u->upper - u->lower);
}

/**
 * Calculates the value of the probability density function for the specified
 * uniform distribution at point x.  If the distribution is NULL then it is
 * assumed to be U[0, 1).
 *
 * @param u  the uniform distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the PDF at point x.
 */
double uniform_pdf(const uniform_distribution u, double x) {
  if (u == NULL)
    return (x < 0.0 || x >= 1.0) ? 0.0 : 1.0;
  return (x < u->lower || x >= u->upper) ? 0.0 : (1.0 / (u->upper - u->lower));
}

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified uniform distribution at point x.  If the
 * distribution is NULL then it is assumed to be U[0, 1).
 *
 * @param u  the uniform distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the first derivative of the PDF at point x.
 */
double uniform_1st_order_pdf(const uniform_distribution u, double x) {
  if (u == NULL)
    return (x < 0.0 || x >= 1.0) ? -HUGE_VAL : 0.0;
  return (x < u->lower || x >= u->upper) ? -HUGE_VAL : 0.0;
}

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified uniform distribution at point x.  If the
 * distribution is NULL then it is assumed to be U[0, 1).
 *
 * @param u  the uniform distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the second derivative of the PDF at point x.
 */
double uniform_2nd_order_pdf(const uniform_distribution u, double x) {
  if (u == NULL)
    return (x < 0.0 || x >= 1.0) ? -HUGE_VAL : 0.0;
  return (x < u->lower || x >= u->upper) ? -HUGE_VAL : 0.0;
}
