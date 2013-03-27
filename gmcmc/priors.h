#ifndef GMCMC_PRIORS_H
#define GMCMC_PRIORS_H

#include <gmcmc/error.h>
#include <gmcmc/rng.h>
#include <stdbool.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Prior type (univariate and continuous).
 */
typedef struct __gmcmc_prior_type_st {
  const char * name;                                    /**< Name of the distribution */
  bool (*init)(void *, va_list);                        /**< Function to initialise and check parameter values */
  double (*sample)(const gmcmc_rng *, const void *);    /**< Sample function */
  double (*evaluate)(double, const void *);             /**< Evaluate function */
  double (*evaluate_1st_order)(double, const void *);   /**< 1st order derivative evaluate function */
  double (*evaluate_2nd_order)(double, const void *);   /**< 2nd order derivative evaluate function */
  size_t n;                                             /**< Amount of memory needed to store parameters */
} gmcmc_prior_type;

/*!
 * Prior (univariate and continuous).
 */
typedef struct __gmcmc_prior_st gmcmc_prior;

/*!
 * Creates a prior distribution.
 *
 * @param [out] p       the prior distribution to create
 * @param [in]  type    the type of prior distribution to create
 * @param [in]  params  distribution parameters
 *
 * @return GMCMC_SUCCESS if the prior was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *           the parameter vector,
 *         GMCMC_ERROR_INVALID_ARGUMENT if the prior type is invalid or if one
 *           of the distribution parameters is invalid.
 */
gmcmc_error gmcmc_prior_create(gmcmc_prior **, const gmcmc_prior_type *, ...);

/*!
 * Creates a prior distribution which is a duplicate of an existing prior.
 *
 * @param [out] p  the prior distribution to create
 * @param [in]  q  the prior distribution to copy
 *
 * @return GMCMC_SUCCESS if the prior was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the parameter vector.
 */
gmcmc_error gmcmc_prior_create_copy(gmcmc_prior **, const gmcmc_prior *);

/*!
 * Copies a prior distribution into an existing prior.  If the source and
 * destination priors have types which have different numbers of parameters this
 * function will resize the destination prior parameter vector.
 * 
 * @param [out] p  the destination prior distribution
 * @param [in]  q  the source prior distribution
 *
 * @return GMCMC_SUCCESS if the prior was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the parameter vector.
 */
gmcmc_error gmcmc_prior_copy(gmcmc_prior *, const gmcmc_prior *);

/*!
 * Destroys a prior distribution.
 *
 * @param [in] p the prior distribution
 */
void gmcmc_prior_destroy(gmcmc_prior *);

/*!
 * Generates a sample from a prior distribution.
 *
 * @param [in] p  a prior distribution
 * @param [in] r  the random number generator to use
 *
 * @return a sample from the prior
 * 
 * @see gmcmc_rng
 */
double gmcmc_prior_sample(const gmcmc_prior *, const gmcmc_rng *);

/*!
 * Evaluates a prior distribution probability density function.
 *
 * @param [in] p  a prior distribution
 * @param [in] x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double gmcmc_prior_evaluate(const gmcmc_prior *, double);

/*!
 * Evaluates the first order derivative of a prior distribution probability
 * density function.
 *
 * @param [in] p  a prior distribution
 * @param [in] x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double gmcmc_prior_evaluate_1st_order(const gmcmc_prior *, double);

/*!
 * Evaluates the second order derivative of a prior distribution probability
 * density function.
 *
 * @param [in] p  a prior distribution
 * @param [in] x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double gmcmc_prior_evaluate_2nd_order(const gmcmc_prior *, double);

// Predefined priors

/*! A uniform prior over (a, b) */
extern const gmcmc_prior_type * gmcmc_prior_uniform;

/*! A gamma prior with parameters (a, b) */
extern const gmcmc_prior_type * gmcmc_prior_gamma;

/*! A normal prior with parameters for mean and standard deviation */
extern const gmcmc_prior_type * gmcmc_prior_normal;

/*! A lognormal prior with parameters for log-mean and standard deviation */
extern const gmcmc_prior_type * gmcmc_prior_lognormal;

#ifdef __cplusplus
}
#endif

#endif
