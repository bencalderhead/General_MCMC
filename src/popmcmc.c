#include <gmcmc/model.h>

/*!
 * Performs Differential Geometric Population MCMC.
 *
 * @param [in]  options  parameter object
 * @param [in]  model    the model
 *
 * @return GMCMC_SUCCESS on success.
 * @see gmcmc_popmcmc_options
 */
gmcmc_error gmcmc_popmcmc(const gmcmc_popmcmc_options * opts, const gmcmc_model * model) {
  // Check options
  if (opts->output == NULL ||
      opts->temperatures == NULL ||
      opts->num_temps <= 0 ||
      islessequal(opts->adapt_rate, 0.0) ||
      islessequal(opts->lower, 0.0) ||
      islessequal(opts->upper, opts->lower)) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return GMCMC_ERROR_INVALID_ARGUMENT;
  }

  // Create a chain for each temperature
  gmcmc_chain ** chains;
  if ((chains = calloc((size_t)opts->num_temps, sizeof(gmcmc_chain *))) == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  const double * params = gmcmc_model_get_params_const(model);
  unsigned int num_params = gmcmc_model_get_num_params(model);
  double stepsize = gmcmc_model_get_stepsize(model);
  unsigned int num_blocks = gmcmc_model_get_num_blocks(model);

  for (int i = 0; i < opts->num_temps; i++) {
    double temperature = gmcmc_model_get_temperature(model, i);
    if ((error = gmcmc_chain_create(&chains[i], num_params, num_blocks, params, temperature, stepsize)) != GMCMC_SUCCESS) {
      for (int j = 0; j < i; j++)
        gmcmc_chain_destroy(chains[j]);
      free(chains);

      GMCMC_ERROR_HANDLER(error);
      return error;
    }
  }

  /*
   * Population MCMC algorithm
   */

  // Main loop
  for (int i = 0; i < opts->burn_in; i++) {
    // Repeat for each chain in the population
    for (int j = 0; j < opts->num_temps; j++) {

      /*
       * MCMC update of parameter values
       */

      // Check if sampling from the prior
      if (gmcmc_chain_get_temperature(chain[j]) == 0.0) {
        // Sample new values from prior

        // Copy old chain values
        gmcmc_chain * old_chain;
        gmcmc_chain_create_copy(&old_chain, chains[j]);

        // Sample new parameter values
        for (int k = 0; k < gmcmc_model_get_num_params(model); k++) {
          gmcmc_prior p;
          gmcmc_model_get_prior(model, k, &p);
          gmcmc_chain_set_param_value(chains[j], k, gmcmc_prior_sample(opts->rng, p));
        }

        // Evaluate the model at new parameters to get LL, gradient, metric, etc.
        gmcmc_error error = gmcmc_model_evaluate(model, chains[k]);

        // Reinstate old values if not successful
        if (error != GMCMC_SUCCESS) {
          gmcmc_chain_copy(chains[j], old_chain);
          gmcmc_chain_destroy(old_chain);
        }
      }
      else {

        gmcmc_block ** blocks;
        if (gmcmc_model_get_blocking_type(model) == GMCMC_BLOCKING_TYPE_RANDOM) {
          // Set up random blocks
          int * idx_temp;
          if ((idx_temp = malloc(gmcmc_model_get_num_params(model) * sizeof(int))) == NULL) {
            error = GMCMC_ERROR_OUT_OF_MEMORY;
            goto end;
          }
        }
      }
    }

    // Exchange the chains between temperatures
    for (int j = 0; j < 3; j++) {
      // Loop through each chain in turn and attempt exchange with chain above
      for (int k = 0; k < opts->num_temps - 1; k++) {
        // Calculate acceptance ratio between two chains
        double newLL = gmcmc_chain_get_LL(chains[k + 1]) * gmcmc_chain_get_temp(chains[k]) +
                       gmcmc_chain_get_LL(chains[k]) * gmcmc_chain_get_temp(chains[k + 1]);

        double oldLL = gmcmc_chain_get_LL(chains[k]) * gmcmc_chain_get_temp(chains[k]) +
                       gmcmc_chain_get_LL(chains[k + 1]) * gmcmc_chain_get_temp(chains[k + 1]);

        double ratio = newLL - oldLL;

        if (ratio > 0.0 || log(gmcmc_rng_get_double(rng)) < min(0.0, ratio))
          gmcmc_chain_swap(chains[k], chains[k + 1]);
      }
    }
  }

  end:

  // Cleanup
  if (chains != NULL) {
    for (int i = 0; i < opts->num_temps; i++)
      gmcmc_chain_destroy(chains[i]);
  }
  free(chains);

  return error;
}
