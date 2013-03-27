#include <gmcmc/chain.h>
#include <gmcmc/block.h>
#include <gmcmc/error.h>
#include <stdlib.h>

struct __gmcmc_chain_st {
  double temperature;           /*!< chain temperature */
  double * params;              /*!< parameter values */
  double * stepsizes;           /*!< parameter step sizes */
  unsigned int num_params;      /*!< number of parameters */
  unsigned int num_blocks;      /*!< number of blocks/stepsizes */

  unsigned int attempted_mutation, accepted_mutation;        // proposal counters
  unsigned int attempted_exchange, accepted_exchange;
/*
  double log_likelihood;
  double * log_prior;
  double * gradient_log_likelihood;
  double * gradient_log_prior;
  double * fisher_info;
  double * hessian_log_prior;

  gmcmc_block * block;          /*!< Current block */*/
};


gmcmc_error gmcmc_chain_create(gmcmc_chain ** chain, unsigned int num_params, unsigned int num_blocks, const double * params, double temperature, double stepsize) {
  if (num_params < num_blocks || params == NULL || isless(temperature, 0.0) || stepsize == 0.0) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return GMCMC_ERROR_INVALID_ARGUMENT;
  }

  if ((*chain = malloc(sizeof(struct __gmcmc_chain_st))) == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  if (((*chain)->params = malloc(num_params * sizeof(double))) == NULL) {
    free(*chain);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY);
  }

  if (((*chain)->stepsizes = malloc(num_blocks * sizeof(double))) == NULL) {
    free((*chain)->params);
    free(*chain);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY);
  }

  for (int i = 0; i < num_params; i++)
    (*chain)->params[i] = params[i];

  for (int i = 0; i < num_blocks; i++)
    (*chain)->stepsizes = stepsize;

  (*chain)->temperature = temperature;
  (*chain)->num_params = num_params;
  (*chain)->num_blocks = num_blocks;
  (*chain)->attempted_mutation = 0;
  (*chain)->accepted_mutation = 0;
  (*chain)->attempted_exchange = 0;
  (*chain)->accepted_exchange = 0;

  return GMCMC_SUCCESS;
}

void gmcmc_chain_destroy(gmcmc_chain * chain) {
  free(chain->params);
  free(chain->stepsizes);
  free(chain);
}

double gmcmc_chain_get_temperature(const gmcmc_chain * chain) {
  return chain->temperature;
}

double gmcmc_chain_get_stepsize(const gmcmc_chain * chain, int i) {
  if (i >= chain->num_blocks) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return 0.0;
  }
  return chain->stepsizes[i];
}

double gmcmc_chain_get_parameter(const gmcmc_chain * chain, int i) {
  if (i >= chain->num_params) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return 0.0;
  }
  return chain->params[i];
}

gmcmc_error gmcmc_chain_set_parameter(gmcmc_chain * chain, int i, double value) {
  if (i >= chain->num_params) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return GMCMC_ERROR_INVALID_ARGUMENT;
  }
  chain->params[i] = value;
  return GMCMC_SUCCESS;
}
