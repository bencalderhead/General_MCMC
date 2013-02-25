#ifndef LOGNORMAL_H
#define LOGNORMAL_H

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
typedef struct __lognormal_distribution_st * lognormal_distribution;

/**
 * Creates a new lognormal distribution.
 *
 * @param l      the newly created lognormal distribution is returned through
 *                 this pointer
 * @param shape  the shape
 * @param scale  the log-scale
 *
 * @return 0 on success, EINVAL if shape <= 0, ENOMEM if there is not enough
 *  memory available to create another lognormal distribution.
 */
int lognormal_create(lognormal_distribution *, double, double);

/**
 * Destroys a lognormal distribution, freeing any memory allocated.
 *
 * @param l  the lognormal distribution to destroy
 */
void lognormal_destroy(lognormal_distribution);

double lognormal_get_shape(const lognormal_distribution);
double lognormal_get_scale(const lognormal_distribution);

int lognormal_set_shape(lognormal_distribution, double);
int lognormal_set_scale(lognormal_distribution, double);

/**
 * Generates a lognormal distributed double precision floating point variable.
 * If the distribution is NULL it is taken to be lnN(1, 1).
 *
 * @param l      the lognormal distribution (may be NULL)
 * @param rng    a function generating double precision floating point variables
 *                 uniformly distributed over [0, 1).
 * @param state  state for the rng function
 *
 * @return a lognormal distributed double precision floating point variable.
 */
double lognormal_rand(const lognormal_distribution, double (*)(const void *), const void *);

/**
 * Calculates the value of the probability density function for the specified
 * lognormal distribution at point x.  If the distribution is NULL it is taken
 * to be lnN(1, 1).
 *
 * @param l  the lognormal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the PDF at point x.
 */
double lognormal_pdf(const lognormal_distribution, double);

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified lognormal distribution at point x.  If the
 * distribution is NULL it is taken to be lnN(1, 1).
 *
 * @param l  the lognormal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the first derivative of the PDF at point x.
 */
double lognormal_1st_order_pdf(const lognormal_distribution, double);

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified lognormal distribution at point x.  If the
 * distribution is NULL it is taken to be lnN(1, 1).
 *
 * @param l  the lognormal distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the second derivative of the PDF at point x.
 */
double lognormal_2nd_order_pdf(const lognormal_distribution, double);

#ifdef __cplusplus
}
#endif

#endif
