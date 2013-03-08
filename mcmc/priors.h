#ifndef PRIORS_H
#define PRIORS_H

#include <stdbool.h>
#include <stdarg.h>
#include "rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Prior type (univariate and continuous).
 */
typedef struct __mcmc_prior_type_st {
  const char * name;                                    // Name of the distribution
  bool (*init)(void *, va_list);                        // Function to initialise and check parameter values
  double (*sample)(const mcmc_rng, const void *);       // Sample function
  double (*evaluate)(double, const void *);             // Evaluate function
  double (*evaluate_1st_order)(double, const void *);   // 1st order derivative evaluate function
  double (*evaluate_2nd_order)(double, const void *);   // 2nd order derivative evaluate function
  size_t n;                                             // Amount of memory needed to store parameters
} * mcmc_prior_type;

/**
 * Prior (univariate and continuous).
 */
typedef struct __mcmc_prior_st * mcmc_prior;

/**
 * Creates a prior distribution.
 *
 * @param p       the prior distribution to create
 * @param type    the type of prior distribution to create
 * @param params  distribution parameters
 *
 * @return MCMC_SUCCESS if the prior was created successfully,
 *         MCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the parameter vector,
 *         MCMC_ERROR_INVALID_ARGUMENT if one of the distribution parameters
 *         is invalid.
 */
mcmc_error mcmc_prior_create(mcmc_prior *, const mcmc_prior_type, ...);

/**
 * Destroys a prior distribution.
 *
 * @param p the prior distribution
 */
void mcmc_prior_destroy(mcmc_prior);

/**
 * Generates a sample from a prior distribution.
 *
 * @param p  a prior distribution
 * @param r  the random number generator to use
 *
 * @return a sample from the prior
 */
double mcmc_prior_sample(const mcmc_prior, const mcmc_rng);

/**
 * Evaluates a prior distribution probability density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double mcmc_prior_evaluate(const mcmc_prior, double);

/**
 * Evaluates the first order derivative of a prior distribution probability
 * density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double mcmc_prior_evaluate_1st_order(const mcmc_prior, double);

/**
 * Evaluates the second order derivative of a prior distribution probability
 * density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double mcmc_prior_evaluate_2nd_order(const mcmc_prior, double);

// Predefined priors

/** A uniform prior over (a, b) */
extern const mcmc_prior_type mcmc_prior_uniform;

/** A gamma prior with parameters (a, b) */
extern const mcmc_prior_type mcmc_prior_gamma;

/** A normal prior with parameters for mean and standard deviation */
extern const mcmc_prior_type mcmc_prior_normal;

/** A lognormal prior with parameters for log-mean and standard deviation */
extern const mcmc_prior_type mcmc_prior_lognormal;

#ifdef __cplusplus
}
#endif

#endif
