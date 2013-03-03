#ifndef MCMC_H
#define MCMC_H

#include "priors.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  FILE * output;                // output file for saving samples
  bool save_burn_in;            // whether to save burn-in
  double * temperatures;        // temperature schedule
  int n_temps;                  // number of tempered distributions to use
  
  int burn_in;                  // number of burn-in samples
  int posterior;                // number of posterior samples
  int save_interval;            // how often to save

  double adapt_rate;            // iteration interval for adapting step sizes
  double upper, lower;          // upper and lower step sizes
  
  int cores;                    // number of cores to use
} mcmc_options;

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
  double * fischer_info;
  double * hessian_log_prior;
  
  int * current_block;
  int current_block_size;
  int current_block_num;
} markov_chain;

typedef enum { BLOCKING_FIXED, BLOCKING_RANDOM } blocking_type;

typedef struct {
  const char * name;

  int (*evaluate)(const model *, markov_chain *);
  int (*proposal)(const model *, markov_chain *, double *, double *);

  const char ** param_names;
  double * params;
  prior * priors;
  double step_size;
  int n;

  struct {
    blocking_type type;
    int * blocks;
    int n;
  } blocking;

  double ** data;
  double * time_points;
  int n_time_points;
} model;

#ifdef __cplusplus
}
#endif

#endif
