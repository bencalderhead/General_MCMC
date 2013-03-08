#include "mcmc/rng.h"

#include <stdlib.h>

/**
 * Random number generator.
 */
struct __mcmc_rng_st {
  mcmc_rng_type type;        // The type of PRNG algorithm used
  void * state;                 // RNG state vector
};

/**
 * Creates a new random number generator
 *
 * @param r     the RNG to create
 * @param type  the RNG algorithm to use
 * @param seed  the seed to initialise the RNG state with
 *
 * @return MCMC_SUCCESS if the RNG was created successfully,
 *         MCMC_ERROR_INVALID_ARGUMENT if an invalid rng type was passed to the
 *                   function,
 *         MCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *                   the RNG or state vector.
 */
mcmc_error mcmc_rng_create(mcmc_rng * r, const mcmc_rng_type type, uint64_t seed) {
  // Check the rng type is valid
  if (type->name == NULL ||
      type->max <= type->min ||
      type->set == NULL ||
      type->get == NULL ||
      type->get_double == NULL)
    return MCMC_ERROR_INVALID_ARGUMENT;

  if ((*r = malloc(sizeof(struct __mcmc_rng_st))) == NULL)
    return MCMC_ERROR_OUT_OF_MEMORY;

  // Allocate space for the state vector
  if (((*r)->state = malloc(type->size)) == NULL) {
    free(*r);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  (*r)->type = type;       // Store a pointer to the type of RNG created
  mcmc_rng_set(*r, seed);  // Seed the RNG to initialise the state vector

  return MCMC_SUCCESS;
}

/**
 * Destroys an RNG.
 *
 * @param r  the RNG to destroy
 */
void mcmc_rng_destroy(mcmc_rng r) {
  free(r->state);       // Free the state vector
  free(r);
}

/**
 * (Re)seeds a random number generator.
 *
 * @param r     the RNG to (re)seed
 * @param seed  the new seed
 */
void mcmc_rng_set(mcmc_rng r, uint64_t seed) {
  r->type->set(r->state, seed);
}

/**
 * Generates a random integer uniformly distributed over the range [min, max].
 *
 * @param r  the random number generator
 * @return a random integer.
 */
uint64_t mcmc_rng_get(const mcmc_rng r) {
  return r->type->get(r->state);
}

/**
 * Generates a random double precision real floating point value uniformly
 * distributed over the range (0, 1).
 *
 * @param r  the random number generator
 * @return a random real.
 */
double mcmc_rng_get_double(const mcmc_rng r) {
  return r->type->get_double(r->state);
}
