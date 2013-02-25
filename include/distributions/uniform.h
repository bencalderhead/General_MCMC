#ifndef UNIFORM_H
#define UNIFORM_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Uniform distribution.
 *
 * Parameters:
 *   lower < upper
 *
 * http://en.wikipedia.org/wiki/Uniform_distribution
 */
typedef struct __uniform_distribution_st * uniform_distribution;

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
int uniform_create(uniform_distribution *, double, double);

/**
 * Destroys a uniform distribution, freeing any memory allocated.
 *
 * @param u  the uniform distribution to destroy
 */
void uniform_destroy(uniform_distribution);

double uniform_get_lower(const uniform_distribution);
double uniform_get_upper(const uniform_distribution);

int uniform_set_lower(uniform_distribution, double);
int uniform_set_upper(uniform_distribution, double);

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
double uniform_rand(const uniform_distribution, double (*)(const void *), const void *);

/**
 * Calculates the value of the probability density function for the specified
 * uniform distribution at point x.  If the distribution is NULL then it is
 * assumed to be U[0, 1).
 *
 * @param u  the uniform distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the PDF at point x.
 */
double uniform_pdf(const uniform_distribution, double);

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified uniform distribution at point x.  If the
 * distribution is NULL then it is assumed to be U[0, 1).
 *
 * @param u  the uniform distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the first derivative of the PDF at point x.
 */
double uniform_1st_order_pdf(const uniform_distribution, double);

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified uniform distribution at point x.  If the
 * distribution is NULL then it is assumed to be U[0, 1).
 *
 * @param u  the uniform distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the second derivative of the PDF at point x.
 */
double uniform_2nd_order_pdf(const uniform_distribution, double);

#ifdef __cplusplus
}
#endif

#endif
