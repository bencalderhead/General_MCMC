#include <gmcmc/parameter.h>
#include <stdlib.h>

/*!
 * Geometric MCMC parameter type.
 */
struct __gmcmc_parameter_st {
  char * name;
  double value;
  gmcmc_prior * prior;
};

/*!
 * Creates a new Geometric MCMC parameter.
 *
 * @param [out] p      the parameter
 * @param [in]  name   a name for the parameter
 * @param [in]  value  the initial value of the parameter
 * @param [in]  prior  the prior distribution for the parameter
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_INVALID_ARGUMENT if the prior is NULL,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create
 *               another parameter.
 */
gmcmc_error gmcmc_parameter_create(gmcmc_parameter ** p, const char * name,
                                   double value, const gmcmc_prior * prior) {
  if (prior == NULL)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);

  if ((*p = malloc(sizeof(gmcmc_parameter))) == NULL)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);

  size_t namelen = strlen(name);
  if (((*p)->name = malloc((namelen + 1) * sizeof(char))) == NULL) {
    free(*p);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
  }

  (*p)->name = strncpy((*p)->name, name, namelen + 1);
  (*p)->value = value;

  gmcmc_error error;
  if ((error = gmcmc_prior_create_copy(&(*p)->prior, prior)) != GMCMC_SUCCESS) {
    free((*p)->name);
    free(*p);
    GMCMC_ERROR_HANDLER(error);
  }

  return GMCMC_SUCCESS;
}

/*!
 * Creates a new Geometric MCMC parameter which is a copy of an existing
 * parameter.
 *
 * @param [out] dest  the parameter to create
 * @param [in]  src   the parameter to copy
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create
 *               another parameter.
 */
gmcmc_error gmcmc_parameter_create_copy(gmcmc_parameter ** dest,
                                        const gmcmc_parameter * src) {
  if ((*dest = malloc(sizeof(gmcmc_parameter))) == NULL)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);

  size_t namelen = strlen(src->name);
  if (((*dest)->name = malloc((namelen + 1) * sizeof(char))) == NULL) {
    free(*dest);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
  }

  (*dest)->name = strncpy((*dest)->name, src->name, namelen + 1);
  (*dest)->value = src->value;

  gmcmc_error error;
  if ((error = gmcmc_prior_create_copy(&(*dest)->prior, src->prior)) !=
    GMCMC_SUCCESS) {
    free((*dest)->name);
    free(*dest);
    GMCMC_ERROR_HANDLER(error);
  }

  return GMCMC_SUCCESS;
}

/*!
 * Copies a Geometric MCMC parameter.
 *
 * @param [out] dest  the destination of the copy
 * @param [in]  src   the parameter to copy
 *
 * @return GMCMC_SUCCESS on success.
 */
gmcmc_error gmcmc_parameter_copy(gmcmc_parameter * dest,
                                 const gmcmc_parameter * src) {
  size_t srcnamelen = strlen(src->name);
  size_t destnamelen = strlen(dest->name);
  if (srcnamelen != destnamelen) {
    char * name;
    if ((name = malloc((srcnamelen + 1) * sizeof(char))) == NULL)
      GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    free(dest->name);
    dest->name = name;
  }

  dest->name = strncpy(dest->name, src->name, srcnamelen + 1);
  dest->value = src->value;

  gmcmc_error error;
  if ((error = gmcmc_prior_copy(&dest->prior, src->prior)) !=
    GMCMC_SUCCESS)
    GMCMC_ERROR_HANDLER(error);

  return GMCMC_SUCCESS;
}

/*!
 * Destroys a Geometric MCMC parameter.
 *
 * @param [in] p  the parameter to destroy
 */
void gmcmc_parameter_destroy(gmcmc_parameter * p) {
  free(p->name);
  free(p);
}

/*!
 * Gets the name of the parameter.
 *
 * @param [in] p  the parameter
 *
 * @return the name.
 */
const char * gmcmc_parameter_get_name(const gmcmc_parameter * p) {
  return p->name;
}

/*!
 * Gets the value associated with a parameter.
 *
 * @param [in] p  the parameter
 *
 * @return the value.
 */
double gmcmc_parameter_get_value(const gmcmc_parameter * p) {
  return p->value;
}

/*!
 * Updates the parameter value with a sample from the prior.
 *
 * @param [in,out] p    the parameter
 * @param [in]     rng  the random number generator to use to generate the
 *                        random sample from the prior
 */
void gmcmc_parameter_update(gmcmc_parameter * p, const gmcmc_rng * rng) {
  p->value = gmcmc_prior_sample(p->prior, rng);
}
