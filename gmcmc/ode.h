#ifndef ODE_H
#define ODE_H

#include "mcmc.h"

typedef struct __ode_model_st * ode_model;

int ode_evaluate_mh(const ode_model, markov_chain);
int ode_evaluate_mh_obs(const ode_model, markov_chain);
int ode_evaluate_smmala(const ode_model, markov_chain);
int ode_evaluate_smmala_numerical(const ode_model, markov_chain);
int ode_evaluate_smmala_numerical_lowrank(const ode_model, markov_chain);
int ode_evaluate_smmala_obs_numerical(const ode_model, markov_chain);

int ode_proposal_mh(const ode_model, const markov_chain, double *, double *);
int ode_proposal_smmala(const ode_model, const markov_chain, double *, double *);

#endif
