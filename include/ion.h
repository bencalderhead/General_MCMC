#ifndef ION_H
#define ION_H

#include "mcmc.h"

/**
 * ION model structure.
 */
struct __ion_model_st {
  model m;              // Every ION Channel model is also a MCMC model

  int n_open, n_closed; // Number of states (open and closed)
  int * open, * closed; // Which species are observed and unobserved (indices into Q matrix)
};

int ion_evaluate_mh(const ion_model *, markov_chain *);
int ion_proposal_mh(const ion_model *, const markov_chain *, double *, double *);

#endif
