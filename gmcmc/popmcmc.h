#ifndef GMCMC_POPMCMC_H
#define GMCMC_POPMCMC_H

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Geometric Population MCMC Options parameter object.
 */
typedef struct {
  gmcmc_rng * rng;              /**< parallel PRNG to use */

  FILE * output;                /**< output file for saving samples */
  bool save_burn_in;            /**< whether to save burn-in */
  double * temperatures;        /**< temperature schedule */
  int num_temps;                /**< number of tempered distributions to use */

  int burn_in;                  /**< number of burn-in samples */
  int posterior;                /**< number of posterior samples */

  double adapt_rate;            /**< iteration interval for adapting step sizes */
  double upper, lower;          /**< upper and lower step sizes */
} gmcmc_popmcmc_options;

/*!
 * Performs Differential Geometric Population MCMC.
 *
 * @param [in]  options  parameter object
 * @param [in]  model    the model
 *
 * @return GMCMC_SUCCESS on success.
 * @see gmcmc_popmcmc_options
 */
gmcmc_error gmcmc_popmcmc(const gmcmc_popmcmc_options *, const gmcmc_model *);

#ifdef __cplusplus
}
#endif

#endif
