#ifndef ION_H
#define ION_H

#include "mcmc.h"

typedef struct __ion_workspace_st ion_workspace;

/**
 * ION model structure.
 */
typedef struct __ion_model_st {
  mcmc_model m;              // Every ION Channel model is also a MCMC model

  unsigned int open, closed;          // Number of states (open and closed)
  void (*update_q)(const double *, double *, size_t);   // Function to update Q matrix based on current parameter values

  ion_workspace * workspace;
} ion_model;

mcmc_error ion_model_create(ion_model * restrict, unsigned int, unsigned int, void (*)(const double *, double *, size_t));
void ion_model_destroy(ion_model *);

/**
 * Calculates the log likelihood of the ion channel data for the current parameter values.
 *
 * @param model   the ion channel model
 * @param params  the current parameter values
 *
 * @return the value of the log likelihood.
 */
mcmc_error ion_evaluate_mh(const ion_model *, markov_chain *);

/**
 * Calculates the proposal mean and covariance based on the current parameter
 * values.
 *
 * @param n         the number of parameters
 * @param params    the current parameter values
 * @param stepsize  the parameter step size
 * @param mean      set to the proposal mean
 * @param cov       set to the proposal covariance
 */
mcmc_error ion_proposal_mh(const ion_model *, const markov_chain *, double *, double *);

#endif
