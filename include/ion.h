#ifndef ION_H
#define ION_H

#include "mcmc.h"

/**
 * ION model structure.
 */
struct __ion_model_st {
  model m;              // Every ION Channel model is also a MCMC model

  int open, closed;     // Number of states (open and closed)
  
  double * Q; size_t ldq;       // Q matrix (open + closed by open + closed)
  double * eq_states;           // equilibrium state vector (open + closed)
  double (*update_q)(const double *, double *, size_t); // Function to update Q matrix based on current parameter values
  
  double * X, * X_inv;          // X and X_inv are max(open,closed) by max(open,closed) matrices
  double * v_Q_FF, * specMat_Q_FF;      // Eigenvalues and spectral matrices of Q_FF
  size_t ldqff;
};

int ion_evaluate_mh(const ion_model *, markov_chain *);
int ion_proposal_mh(const ion_model *, const markov_chain *, double *, double *);

#endif
