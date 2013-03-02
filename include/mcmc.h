#ifndef MCMC_H
#define MCMC_H

/**
 * A probability distribution.
 */
typedef struct {
  /**
   * A function to generate a sample from the distribution.
   *
   * @param r       a random number generator to use to generate the sample
   * @param params  distribution parameters
   *
   * @return a sample from the distribution, given the parameters
   */
  double (*rand)(const rng * restrict , const void * restrict);

  /**
   * Probability density function.
   *
   * @param x       the point at which to evaluate the function
   * @param params  distribution parameters
   *
   * @return the value of the probability density function at x.
   */
  double (*pdf)(double, const void *);

  /**
   * Distribution parameters
   */
  void * params;
} prior;

static inline double prior_rand();

typedef enum { BLOCKING_FIXED, BLOCKING_RANDOM } mcmc_blocking_type;

/**
 * MCMC model structure.
 */
struct __model_st {
  const char * type;    // Model type
  const char * name;    // Model name
  const char * sampler; // Model sampler

  int  (*evaluate)(const model *, markov_chain *, int *);       // Evaluation function
  void (*proposal)(const model *, const markov_chain *, double *, double *);    // Proposal function

  const char ** param_names;    // Parameter names
  int n_params;                 // Number of parameters
  bool log10_space;             // Whether to use log space for parameter values
  double * params;              // Starting values
  double step_size;             // Step size

  mcmc_blocking_type blocking;  // Either 'fixed' or 'random'
  int n_blocks;                 // Number of blocks
  int * blocks;                 // Must sum to number of parameters

  prior * priors;               // Prior distribution for each parameter

  double * data;                // Model data
  double * time_points;         // Model timepoints
  int n_time_points;            // Number of timepoints
};
typedef struct __markov_chain_st markov_chain;

#endif
