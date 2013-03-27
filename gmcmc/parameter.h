#ifndef GMCMC_PARAMETER_H
#define GMCMC_PARAMETER_H

#include <gmcmc/priors.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Geometric MCMC parameter type.
 */
typedef struct __gmcmc_parameter_st gmcmc_parameter;

/*!
 * Creates a new Geometric MCMC parameter.
 *
 * @param [out] p      the parameter
 * @param [in]  name   a name for the parameter
 * @param [in]  value  the initial value of the parameter
 * @param [in]  prior  the prior distribution for the parameter
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create
 *               another parameter.
 */
gmcmc_error gmcmc_parameter_create(gmcmc_parameter **, const char *, double, const gmcmc_prior *);

/*!
 * Creates a new Geometric MCMC parameter which is a copy of an existing parameter.
 *
 * @param [out] dest  the parameter to create
 * @param [in]  src   the parameter to copy
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create
 *               another parameter.
 */
gmcmc_error gmcmc_parameter_create_copy(gmcmc_parameter **, const gmcmc_parameter *);

/*!
 * Copies a Geometric MCMC parameter.
 *
 * @param [out] dest  the destination of the copy
 * @param [in]  src   the parameter to copy
 *
 * @return GMCMC_SUCCESS on success.
 */
gmcmc_error gmcmc_parameter_copy(gmcmc_parameter *, const gmcmc_parameter *);

/*!
 * Destroys a Geometric MCMC parameter.
 *
 * @param [in] p  the parameter to destroy
 */
void gmcmc_parameter_destroy(gmcmc_parameter *);

/*!
 * Gets the name of the parameter.
 *
 * @param [in] p  the parameter
 *
 * @return the name.
 */
const char * gmcmc_parameter_get_name(const gmcmc_parameter *);

/*!
 * Gets the value associated with a parameter.
 *
 * @param [in] p  the parameter
 *
 * @return the value.
 */
double gmcmc_parameter_get_value(const gmcmc_parameter *);

/*!
 * Updates the parameter value with a sample from the prior.
 *
 * @param [in,out] p  the parameter
 */
void gmcmc_parameter_update(gmcmc_parameter *);

#ifndef __cplusplus
}
#endif

#endif
