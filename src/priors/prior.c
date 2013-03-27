#include <gmcmc/priors.h>
#include <stdlib.h>
#include <string.h>

/**
 * Prior (univariate and continuous).
 */
struct __gmcmc_prior_st {
  gmcmc_prior_type * type;      /**< Type of prior */
  void * params;                /**< Parameter vector */
};

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
gmcmc_error gmcmc_prior_create(gmcmc_prior ** p, const gmcmc_prior_type * type, ...) {
  // Check for a valid prior type
  if (type->name == NULL || type->init == NULL || type->sample == NULL ||
      type->evaluate == NULL || type->evaluate_1st_order == NULL ||
      type->evaluate_2nd_order == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return GMCMC_ERROR_INVALID_ARGUMENT;
  }

  // Allocate space for the prior
  if ((*p = calloc(1, sizeof(struct __gmcmc_prior_st))) == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate space for a copy of the prior type and for the parameter vector
  if (((*p)->type = malloc(sizeof(gmcmc_prior_type))) == NULL ||
      ((*p)->params = malloc(type->n)) == NULL) {
    free((*p)->type);
    free((*p)->params);
    free(*p);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Store a pointer to the type of prior created
  (*p)->type = memcpy((*p)->type, type, sizeof(gmcmc_prior_type));

  // Initialise the parameter vector
  va_list list;
  va_start(list, type);
  bool valid = type->init((*p)->params, list);
  va_end(list);

  // Check parameter values are valid
  if (!valid) {
    free((*p)->type);
    free((*p)->params);
    free(*p);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return GMCMC_ERROR_INVALID_ARGUMENT;
  }

  return GMCMC_SUCCESS;
}

/*!
 * Creates a prior distribution which is a duplicate of an existing prior.
 *
 * @param [out] p    the prior distribution to create
 * @param [in]  src  the prior distribution to copy
 *
 * @return GMCMC_SUCCESS if the prior was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the parameter vector.
 */
gmcmc_error gmcmc_prior_create_copy(gmcmc_prior ** p, const gmcmc_prior * src) {
  // Allocate space for the prior
  if ((*p = calloc(1, sizeof(struct __gmcmc_prior_st))) == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate space for a copy of the type and parameter vector
  if (((*p)->type = malloc(sizeof(gmcmc_prior_type))) == NULL ||
      ((*p)->params = malloc(src->type->n)) == NULL) {
    free((*p)->type);
    free((*p)->params);
    free(*p);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Copy the type and parameter vector
  (*p)->type = memcpy((*p)->type, src->type, sizeof(gmcmc_prior_type));
  (*p)->params = memcpy((*p)->params, src->params, src->type->n);

  return GMCMC_SUCCESS;
}

/*!
 * Copies a prior distribution into an existing prior.  If the source and
 * destination priors have types which have different numbers of parameters this
 * function will resize the destination prior parameter vector.
 * 
 * @param [out] p    the destination prior distribution
 * @param [in]  src  the source prior distribution
 *
 * @return GMCMC_SUCCESS if the prior was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the parameter vector.
 */
gmcmc_error gmcmc_prior_copy(gmcmc_prior * p, const gmcmc_prior * src) {
  // If source == destination return now
  if (p == src)
    return GMCMC_SUCCESS;
  
  // If the priors have different sizes of parameter vectors
  if (p->type->n != src->type->n) {
    // Allocate space for new parameter vector
    void * params;
    if ((params = malloc(src->type->n)) == NULL) {
      GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
      return GMCMC_ERROR_OUT_OF_MEMORY;
    }

    free(p->params);    // Free existing parameter vector
    p->params = params; // Update parameter vector
  }

  // Copy the type and parameter vector
  p->type = memcpy(p->type, src->type, sizeof(gmcmc_prior_type));
  p->params = memcpy(p->params, src->params, src->type->n);

  return GMCMC_SUCCESS;
}

/*!
 * Destroys a prior distribution.
 *
 * @param [in] p the prior distribution
 */
void gmcmc_prior_destroy(gmcmc_prior * p) {
  free(p->type);
  free(p->params);
  free(p);
}

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
double gmcmc_prior_sample(const gmcmc_prior * p, const gmcmc_rng * r) {
  return p->type->sample(r, p->params);
}

/*!
 * Evaluates a prior distribution probability density function.
 *
 * @param [in] p  a prior distribution
 * @param [in] x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double gmcmc_prior_evaluate(const gmcmc_prior * p, double x) {
  return p->type->evaluate(x, p->params);
}

/*!
 * Evaluates the first order derivative of a prior distribution probability
 * density function.
 *
 * @param [in] p  a prior distribution
 * @param [in] x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double gmcmc_prior_evaluate_1st_order(const gmcmc_prior * p, double x) {
  return p->type->evaluate_1st_order(x, p->params);
}

/*!
 * Evaluates the second order derivative of a prior distribution probability
 * density function.
 *
 * @param [in] p  a prior distribution
 * @param [in] x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
double gmcmc_prior_evaluate_2nd_order(const gmcmc_prior * p, double x) {
  return p->type->evaluate_2nd_order(x, p->params);
}
