#ifndef ION_H
#define ION_H

#include "mcmc.h"

typedef struct __ion_workspace_st ion_workspace;

/**
 * ION model structure.
 */
typedef struct __ion_model_st {
  mcmc_model m;              // Every ION Channel model is also a MCMC model

  int open, closed;          // Number of states (open and closed)
  void (*update_q)(const double *, double *, size_t);   // Function to update Q matrix based on current parameter values

  void * workspace;
} ion_model;

int ion_model_create(ion_model * restrict, int, int, void (*)(const double *, double *, size_t));
void ion_model_destroy(ion_model *);

int ion_evaluate_mh(const ion_model *, markov_chain *);
int ion_proposal_mh(const ion_model *, const markov_chain *, double *, double *);

#endif
