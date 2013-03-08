#include "mcmc/priors.h"

#include <stdlib.h>

/**
 * Prior (univariate and continuous).
 */
struct __mcmc_prior_st {
  mcmc_prior_type type; // Type of prior
  void * params;        // Parameter vector
};

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
mcmc_error mcmc_prior_create(mcmc_prior * p, const mcmc_prior_type type, ...) {
  if ((*p = malloc(sizeof(struct __mcmc_prior_st))) == NULL)
    return MCMC_ERROR_OUT_OF_MEMORY;

  // Allocate space for the parameter vector
  if (((*p)->params = malloc(type->n)) == NULL) {
    free(*p);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  (*p)->type = type;    // Store a pointer to the type of prior created

  // Initialise the parameter vector
  va_list list;
  va_start(list, type);
  bool valid = type->init((*p)->params, list);
  va_end(list);

  // Check parameter values are valid
  if (!valid) {
    free((*p)->params);
    free(*p);
    return MCMC_ERROR_INVALID_ARGUMENT;
  }

  return MCMC_SUCCESS;
}

/**
 * Destroys a prior distribution.
 *
 * @param p the prior distribution
 */
void mcmc_prior_destroy(mcmc_prior p) {
  free(p->params);      // Free the parameter vector
  free(p);
}

/**
 * Generates a sample from a prior distribution.
 *
 * @param p  a prior distribution
 * @param r  the random number generator to use
 *
 * @return a sample from the prior
 */
double mcmc_prior_sample(const mcmc_prior p, const mcmc_rng r) {
  return p->type->sample(r, p->params);
}

/**
 * Evaluates a prior distribution probability density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double mcmc_prior_evaluate(const mcmc_prior p, double x) {
  return p->type->evaluate(x, p->params);
}

/**
 * Evaluates the first order derivative of a prior distribution probability
 * density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double mcmc_prior_evaluate_1st_order(const mcmc_prior p, double x) {
  return p->type->evaluate_1st_order(x, p->params);
}

/**
 * Evaluates the second order derivative of a prior distribution probability
 * density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double mcmc_prior_evaluate_2nd_order(const mcmc_prior p, double x) {
  return p->type->evaluate_2nd_order(x, p->params);
}
