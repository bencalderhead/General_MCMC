#include "mcmc.h"

#include <stdlib.h>

struct __mcmc_model_st {
  const char * name;

  mcmc_evaluate evaluate;
  mcmc_proposal proposal;

  const char ** param_names;
  double * params;
  mcmc_prior * priors;
  int num_params;
  double stepsize;

  mcmc_blocking blocking;

  void * modelspecific;

  double * data;
  double * timepoints;
  int num_timepoints;
};

/**
 * Creates an MCMC model.
 *
 * @param model     the model to create
 * @param name      the name of the model
 * @param evaluate  the model evaluation function
 * @param proposal  the model proposal function
 *
 * @return MCMC_SUCCESS on success,
 *         MCMC_ERROR_INVALID_ARGUMENT if the name, evaluation function or
 *                    proposal function is NULL,
 *         MCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create the
 *                    model.
 */
mcmc_error mcmc_model_create(mcmc_model * model, const char * name,
                             mcmc_evaluate evalulate, mcmc_proposal proposal) {
  if (name == NULL || evaluate == NULL || proposal == NULL)
    return MCMC_ERROR_INVALID_ARGUMENT;

  if ((*model = malloc(sizeof(struct __mcmc_model_st))) == NULL)
    return MCMC_ERROR_OUT_OF_MEMORY;

  size_t length = strlen(name) + 1;
  if (((*model)->name = malloc(length)) == NULL) {
    free(*model);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  (*model)->name = strcpy((*model)->name, name);
  (*model)->evaluate = evaluate;
  (*model)->proposal = proposal;

  // Initialise the rest of the fields to default values
  (*model)->param_names = NULL;
  (*model)->params = NULL;
  (*model)->priors = NULL;
  (*model)->num_params = 0;
  (*model)->stepsize = 1.0;
  (*model)->blocking = NULL;
  (*model)->modelspecific = NULL;
  (*model)->data = NULL;
  (*model)->timepoints = NULL;
  (*model)->num_timepoints = 0;

  return MCMC_SUCCESS;
}

/**
 * Destroys an MCMC model.
 *
 * @param model  the model to destroy
 */
void mcmc_model_destroy(mcmc_model model) {
  free(model->name);
  free(name);
}

// Parameters

/**
 * Sets the parameters in the model.
 *
 * @param model   the model
 * @param names   the names of the parameters (may contain NULL values or be NULL)
 * @param params  the initial values of the parameters (may be NULL)
 * @param priors  prior distributions for each parameter
 * @param n       the number of parameters
 *
 * @return MCMC_SUCCESS on success,
 *         MCMC_ERROR_INVALID_ARGUMENT if n < 0 or n > 0 and priors is NULL
 *         MCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to copy the
 *                    parameters into the model.
 */
mcmc_error mcmc_model_set_parameters(mcmc_model model, const char ** names,
                                     const double * params,
                                     const mcmc_prior * priors, int n) {
  if (n < 0 || (n > 0 && priors == NULL))
    return MCMC_ERROR_INVALID_ARGUMENT;

  // Allocate memory for and create copies of the new parameter values
  const char ** param_names;
  if (names != NULL && n > 0) {
    if ((param_names = malloc((size_t)n * sizeof(const char *))) == NULL)
      return MCMC_ERROR_OUT_OF_MEMORY;

    for (int i = 0; i < n; i++) {
      if (names[i] != NULL) {
        size_t length = strlen(names[i]) + 1;
        if ((param_names[i] = malloc(length)) == NULL) {
          for (int j = i - 1; j >= 0; j--)
            free(param_names[i]);
          free(param_names);
          return MCMC_ERROR_OUT_OF_MEMORY;
        }
        param_names[i] = strcpy(param_names[i], names[i]);
      }
      else
        names[i] = NULL;
    }
  }
  else
    param_names = NULL;

  double * new_params;
  if ((new_params = malloc((size_t)n * sizeof(double))) == NULL) {
    for (int i = 0; i < n; i++)
      free(param_names[i]);
    free(param_names);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }
  new_params = memcpy(new_params, params, (size_t)n * sizeof(double));

  mcmc_prior * new_priors;
  if ((new_priors = malloc((size_t)n * sizeof(mcmc_prior))) == NULL) {
    free(new_params);
    for (int i = 0; i < n; i++)
      free(param_names[i]);
    free(param_names);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }
  for (int i = 0; i < n: i++) {
    mcmc_error error;
    if ((error = mcmc_prior_create_copy(&new_priors[i], priors[i])) != MCMC_SUCCESS) {
      for(int j = i - 1; j >= 0; j--)
        mcmc_prior_destroy(new_priors[i]);
      free(new_priors);
      free(new_params);
      for (int i = 0; i < n; i++)
        free(param_names[i]);
      free(param_names);
      return error;
    }
  }

  // free existing parameters
  if (model->param_names != NULL) {
    for (int i = 0; i < model->n; i++)
      free(model->param_names[i]);
  }
  for (int i = 0; i < model->n; i++)
    mcmc_prior_destroy(model->priors[i]);
  free(model->param_names);
  free(model->params);
  free(model->priors);

  // Replace the parameters with copies of the new ones
  model->param_names = param_names;
  model->params = new_params;
  model->priors = new_priors;
  model->num_params = n;

  return MCMC_SUCCESS;
}

mcmc_error mcmc_model_set_stepsize(mcmc_model model, double stepsize) {
  if (isnan(stepsize) || isinf(stepsize) || islessequal(stepsize, 0.0))
    return MCMC_ERROR_INVALID_ARGUMENT;
  model->stepsize = stepsize;
  return MCMC_SUCCESS;
}

mcmc_error mcmc_model_set_blocking(mcmc_model model, const mcmc_blocking blocking) {
  // Check that for random blocking the block indices sum to the number of parameters
  mcmc_blocking_type type = mcmc_blocking_get_type(blocking);
  if (type == MCMC_BLOCKING_RANDOM) {
    int sum = 0;
    int n = mcmc_blocking_get_num_blocks(blocking);
    for (int i = 0; i < n; i++) {
      const mcmc_block block = mcmc_blocking_get_block(blocking, i);
      int size = mcmc_block_get_size(block);
      for (int j = 0; j < size; j++)
        sum += mcmc_block_get_index(block, j);
    }
    if (size != model->num_params)
      return MCMC_ERROR_INVALID_ARGUMENT;
  }
}

mcmc_error mcmc_model_set_modelspecific(mcmc_model, void *);
mcmc_error mcmc_model_set_data(mcmc_model, const double *, const double *, int);

const char * mcmc_model_get_name(const mcmc_model);
mcmc_error mcmc_model_evaluate(const mcmc_model, markov_chain);
mcmc_error mcmc_model_proposal(const mcmc_model, markov_chain, double *, double *);

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
