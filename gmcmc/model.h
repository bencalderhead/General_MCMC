#ifndef GMCMC_MODEL_H
#define GMCMC_MODEL_H

#include <gmcmc/priors.h>
#include <gmcmc/block.h>
#include <gmcmc/chain.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Geometric MCMC likelihood evaluation function type.
 *
 * @param [in]  model  MCMC model
 * @param [out] chain  Current value in the chain
 *
 * @return error code.
 */
typedef gmcmc_error (*gmcmc_evaluate)(const gmcmc_model *, gmcmc_chain *);

/*!
 * Geometric MCMC proposal function type.
 *
 * @param [in]  model  MCMC model
 * @param [in]  chain  Current value in the chain
 * @param [out] mean   the mean value of the proposal
 * @param [out] var    the variance of the proposal
 *
 * @return error code.
 */
typedef gmcmc_error (*gmcmc_proposal)(const gmcmc_model *, const gmcmc_chain *,
                                      double *, double *);

/*!
 * Geometric MCMC model.
 */
typedef struct __gmcmc_model_st gmcmc_model;

gmcmc_error gmcmc_model_create(gmcmc_model **, gmcmc_evaluate, gmcmc_proposal, const char *);
gmcmc_error gmcmc_model_create_copy(gmcmc_model **, const gmcmc_model *);
gmcmc_error gmmcmc_model_copy(gmcmc_model *, const gmcmc_model *);
void gmcmc_model_destroy(gmcmc_model *);

/*!
 * Set the parameters in the model.
 *
 * @param [out] model   the model
 * @param [in]  name    the names of the parameters (may be NULL or contain NULL)
 * @param [in]  params  initial parameter values (may be NULL)
 * @param [in]  priors  prior distributions for each parameter
 * @param [in]  n       number of parameters
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_INVALID_ARGUMENT if priors is NULL or n is less than zero,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to copy the
 *               parameters into the model.
 */
gmcmc_error gmcmc_model_set_parameters(gmcmc_model *, const char **, const double *, const gmcmc_prior *, int);

/*!
 * Sets the parameter step size.
 *
 * @param [out] model     the model
 * @param [in]  stepsize  the step size
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_INVALID_ARGUMENT if the step size is less than or equal
 *               to zero.
 */
gmcmc_error gmcmc_model_set_stepsize(gmcmc_model *, double);
gmcmc_error gmcmc_model_set_blocking(gmcmc_model *, const gmcmc_block *, int);
gmcmc_error gmcmc_model_set_data(gmcmc_model *, const double *, const double *, int);

const char * gmcmc_model_get_name(const gmcmc_model *);
gmcmc_error gmcmc_model_evaluate(const gmcmc_model *, gmcmc_chain *);
gmcmc_error gmcmc_model_proposal(const gmcmc_model *, gmcmc_chain *, double *, double *);

int mcmc_model_get_num_parameters(const mcmc_model);
const char * mcmc_model_get_parameter_name(const mcmc_model, int);
double mcmc_model_get_parameter_value(const mcmc_model, int);
const mcmc_prior mcmc_model_get_parameter_prior(const mcmc_model, int);
double mcmc_model_get_stepsize(const mcmc_model);
const mcmc_blocking mcmc_model_get_blocking(const mcmc_model);

void * mcmc_model_get_model_specific(const mcmc_model);

const double * mcmc_model_get_data(const mcmc_model);
const double * mcmc_model_get_timepoints(const mcmc_model);
int mcmc_model_get_num_timepoints(const mcmc_model);

#ifdef __cplusplus
}
#endif

#endif
