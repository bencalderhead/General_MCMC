#ifndef MCMC_H
#define MCMC_H

#include "priors.h"
#include "blocking.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  FILE * output;                // output file for saving samples
  bool save_burn_in;            // whether to save burn-in
  double * temperatures;        // temperature schedule
  unsigned int num_temps;       // number of tempered distributions to use

  unsigned int burn_in;         // number of burn-in samples
  unsigned int posterior;       // number of posterior samples
  unsigned int save_interval;   // how often to save

  double adapt_rate;            // iteration interval for adapting step sizes
  double upper, lower;          // upper and lower step sizes

  unsigned int cores;           // number of cores to use
} mcmc_options;


typedef mcmc_error (*mcmc_evaluate)(const mcmc_model, markov_chain);
typedef mcmc_error (*mcmc_proposal)(const mcmc_model, markov_chain, double *, double *);

typedef struct __mcmc_model_st {
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
} * mcmc_model;

mcmc_error mcmc_model_create(mcmc_model *, const char *, mcmc_evaluate, mcmc_proposal);
void mcmc_model_destroy(mcmc_model);

mcmc_error mcmc_model_set_parameters(mcmc_model, const char **, const double *, const mcmc_prior *, int);
mcmc_error mcmc_model_set_stepsize(mcmc_model, double);
mcmc_error mcmc_model_set_blocking(mcmc_model, const mcmc_blocking);
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

typedef struct {
  double temperature;   // chain temperature
  double * params;      // parameter values
  double * stepsize;    // parameter step sizes for parameters in each chain in each population

  int * attempted_mutation, * accepted_mutation;        // proposal counters
  int attempted_exchange, accepted_exchange;

  double * log_likelihood;
  double * log_prior;
  double * gradient_log_likelihood;
  double * gradient_log_prior;
  double * fisher_info;
  double * hessian_log_prior;

  int * current_block;
  int current_block_size;
  int current_block_num;
} markov_chain;

#ifdef __cplusplus
}
#endif

#endif
