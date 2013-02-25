#ifndef GAMMA_H
#define GAMMA_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Gamma distribution.
 *
 * Parameters:
 *   shape > 0
 *   scale > 0
 *
 * Ref: http://en.wikipedia.org/wiki/Gamma_distribution
 */
typedef struct __gamma_distribution_st * gamma_distribution;

/**
 * Creates a new gamma distribution.
 *
 * @param g      the newly created gamma distribution is returned through
 *                 this pointer
 * @param shape  the shape
 * @param scale  the scale
 *
 * @return 0 on success, EINVAL if shape <= 0 or scale <= 0, ENOMEM if there is
 *  not enough memory available to create another gamma distribution.
 */
int gamma_create(gamma_distribution *, double, double);

/**
 * Destroys a gamma distribution, freeing any memory allocated.
 *
 * @param g  the gamma distribution to destroy
 */
void gamma_destroy(gamma_distribution);

double gamma_get_shape(const gamma_distribution);
double gamma_get_scale(const gamma_distribution);

int gamma_set_shape(gamma_distribution, double);
int gamma_set_scale(gamma_distribution, double);

/**
 * Generates a gamma distributed double precision floating point variable.
 * If the distribution is NULL it is taken to be G(1, 1).
 *
 * @param g      the gamma distribution (may be NULL)
 * @param rng    a function generating double precision floating point variables
 *                 uniformly distributed over [0, 1).
 * @param state  state for the rng function
 *
 * @return a gamma distributed double precision floating point variable.
 */
double gamma_rand(const gamma_distribution, double (*)(const void *), const void *);

/**
 * Calculates the value of the probability density function for the specified
 * gamma distribution at point x.  If the distribution is NULL it is taken
 * to be G(1, 1).
 *
 * @param g  the gamma distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the PDF at point x.
 */
double gamma_pdf(double, const gamma_distribution);

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified gamma distribution at point x.  If the
 * distribution is NULL it is taken to be G(1, 1).
 *
 * @param l  the gamma distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the first derivative of the PDF at point x.
 */
double gamma_1st_order_pdf(double, const gamma_distribution);

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified gamma distribution at point x.  If the
 * distribution is NULL it is taken to be G(1, 1).
 *
 * @param l  the gamma distribution (may be NULL)
 * @param x  the point to evaluate the PDF at
 * @return the value of the second derivative of the PDF at point x.
 */
double gamma_2nd_order_pdf(double, const gamma_distribution);

#ifdef __cplusplus
}
#endif

#endif
