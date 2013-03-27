#include <gmcmc/model.h>

struct __gmcmc_model_st {
  const char * name;            /*!< Model name */

  gmcmc_proposal proposal;      /*!< Proposal function */
  gmcmc_evaluate evaluate;      /*!< Evaluate function */

  gmcmc_parameter ** params;    /*!< Parameters */
  int num_params;               /*!< Number of parameters */
  double stepsize;              /*!< Parameter step size */

  gmcmc_blocking_type type;     /*!< Parameter blocking type */
  gmcmc_block * blocks;         /*!< Parameter blocks */
  int num_blocks;               /*!< Number of parameter blocks */

  void * modelspecific;         /*!< Model specific data */

  double * data;                /*!< Data */
  double * timepoints;          /*!< Timepoints */
  int num_timepoints;           /*!< Number of timepoints */
};


const char * gmcmc_model_get_name(const gmcmc_model * model) {
  return model->name;
}

gmcmc_error gmcmc_model_proposal(const gmcmc_model * model, gmcmc_chain * chain, double * mean, double * cov) {
  return model->proposal(model, chain, mean, cov);
}

gmcmc_error gmcmc_model_evaluate(const gmcmc_model * model, gmcmc_chain * chain) {
  return model->evaluate(model, chain);
}

const gmcmc_parameter * gmcmc_model_get_parameter(const gmcmc_model * model, int i) {
  if (i >= model->num_params)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
  return model->parameters[i];
}

gmcmc_error gmcmc_model_set_parameters(gmcmc_model * model, const gmcmc_parameter ** params, int n) {

}

gmcmc_error gmcmc_model_set_stepsize(gmcmc_model * model, double stepsize) {
}

