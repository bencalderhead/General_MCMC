#ifndef RNG_H
#define RNG_H

#include <stddef.h>
#include <stdint.h>

#include "error.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Random number generator algorithm.
 */
typedef struct __mcmc_rng_type_st {
  const char * name;    // Name of the algorithm
  uint64_t min, max;    // Minimum and maximum integer values generated
  size_t size;          // Size of state vector
  void (*set)(void *, uint64_t);  // Seed function
  uint64_t (*get)(void *);        // Integer generation function
  double (*get_double)(void *);   // Real generation function on (0,1)
} * mcmc_rng_type;

/**
 * Random number generator.
 */
typedef struct __mcmc_rng_st * mcmc_rng;

/**
 * Creates a new random number generator
 *
 * @param r     the RNG to create
 * @param type  the RNG algorithm to use
 * @param seed  the seed to initialise the RNG state with
 *
 * @return MCMC_SUCCESS if the RNG was created successfully,
 *         MCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the RNG state vector.
 */
mcmc_error mcmc_rng_create(mcmc_rng *, const mcmc_rng_type, uint64_t);

/**
 * Destroys an RNG.
 *
 * @param r  the RNG to destroy
 */
void mcmc_rng_destroy(mcmc_rng);

/**
 * (Re)seeds a random number generator.
 *
 * @param r     the RNG to (re)seed
 * @param seed  the new seed
 */
void mcmc_rng_set(mcmc_rng, uint64_t);

/**
 * Generates a random integer uniformly distributed over the range [min, max].
 *
 * @param r  the random number generator
 * @return a random integer.
 */
uint64_t mcmc_rng_get(const mcmc_rng);

/**
 * Generates a random double precision real floating point value uniformly
 * distributed over the range (0, 1).
 *
 * @param r  the random number generator
 * @return a random real.
 */
double mcmc_rng_get_double(const mcmc_rng);

// Random number generator algorithms

/** Mersenne Twister 64-bit */
extern const mcmc_rng_type mcmc_rng_mt19337_64;

#ifdef __cplusplus
}
#endif

#endif
