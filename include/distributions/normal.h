#ifndef NORMAL_H
#define NORMAL_H

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
typedef struct __normal_distribution_st * normal_distribution;

/**
 * Creates a new normal distribution.
 *
 * @param n     the newly created normal distribution is returned through this
 *                pointer
 * @param mean  the mean
 * @param sd    the standard deviation
 *
 * @return 0 on success, EINVAL if sd <= 0, ENOMEM if there is not enough
 *  memory available to create another normal distribution.
 */
int normal_create(normal_distribution *, double, double);

/**
 * Destroys a normal distribution, freeing any memory allocated.
 *
 * @param n  the normal distribution to destroy
 */
void normal_destroy(normal_distribution);

double normal_get_mean(const normal_distribution);
double normal_get_variance(const normal_distribution);

int normal_set_mean(normal_distribution, double);
int normal_set_standard_deviation(normal_distribution, double);

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
double normal_rand(const normal_distribution, double (*)(const void *), const void *);

/**
 * Calculates the log of the value of the probability density function for the
 * specified normal distribution at point x.  If the distribution is NULL then
 * it is assumed to be N(0, 1).
 *
 * @param n  the normal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the log of the value of the PDF at point x.
 */
double normal_log_pdf(const normal_distribution, double);

/**
 * Calculates the value of the probability density function for the specified
 * normal distribution at point x.  If the distribution is NULL then it is
 * assumed to be N(0, 1).
 *
 * @param n  the normal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the PDF at point x.
 */
double normal_pdf(const normal_distribution, double);

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified normal distribution at point x.  If the
 * distribution is NULL then it is assumed to be N(0, 1).
 *
 * @param n  the normal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the first derivative of the PDF at point x.
 */
double normal_1st_order_pdf(const normal_distribution, double);

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified normal distribution at point x.  If the
 * distribution is NULL then it is assumed to be N(0, 1).
 *
 * @param n  the normal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the second derivative of the PDF at point x.
 */
double normal_2nd_order_pdf(const normal_distribution, double);

#ifdef __cplusplus
}
#endif

#endif
