#include "gmcmc/rng.h"

#include <stdlib.h>
#include <string.h>

/**
 * Random number generator.
 */
struct __gmcmc_rng_st {
  gmcmc_rng_type * type;        /**< The type of PRNG algorithm used */
  void * state;                 /**< RNG state vector */
};

/*!
 * Creates a new random number generator and initialises it with a seed.
 *
 * @param [out] r     the RNG to create
 * @param [in]  type  the RNG algorithm to use
 * @param [in]  seed  the seed to initialise the RNG state with
 *
 * @return GMCMC_SUCCESS if the RNG was created successfully,
 *         GMCMC_ERROR_INVALID_ARGUMENT if type is not a valid RNG type,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the RNG state vector.
 */
gmcmc_error gmcmc_rng_create(gmcmc_rng ** r, const gmcmc_rng_type * type,
                             uint64_t seed) {
  // Check the RNG type passed
  if (type->name == NULL || type->min >= type->max || type->set == NULL ||
      type->get == NULL || type->get_double == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
    return GMCMC_ERROR_INVALID_ARGUMENT;
  }

  // Allocate space for the RNG
  if ((*r = malloc(sizeof(struct __gmcmc_rng_st))) == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate space for a copy of the type and for the state vector
  if (((*r)->type == malloc(sizeof(gmcmc_rng_type))) == NULL || 
      ((*r)->state = malloc(type->size)) == NULL) {
    free((*r)->type);
    free((*r)->state);
    free(*r);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Copy the type
  (*r)->type = memcpy((*r)->type, type, sizeof(gmcmc_rng_type));
  
  // Seed the RNG to initialise the state vector
  gmcmc_rng_set(*r, seed);

  return GMCMC_SUCCESS;
}

/*!
 * Creates a random number generator which is a duplicate of an existing random
 * number generator.
 *
 * @param [out] r  the random number generator to create
 * @param [in]  s  the random number generator to copy
 *
 * @return GMCMC_SUCCESS if the prior was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the RNG state vector.
 */
gmcmc_error gmcmc_rng_create_copy(gmcmc_rng ** r, const gmcmc_rng * s) {
  if ((*r = malloc(sizeof(struct __gmcmc_rng_st))) == NULL) {
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Allocate space for a copy of the type and for the state vector
  if (((*r)->type == malloc(sizeof(gmcmc_rng_type))) == NULL || 
      ((*r)->state = malloc(s->type->size)) == NULL) {
    free((*r)->type);
    free((*r)->state);
    free(*r);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
    return GMCMC_ERROR_OUT_OF_MEMORY;
  }

  // Copy the type and state vector
  (*r)->type = memcpy((*r)->type, s->type, sizeof(gmcmc_rng_type));
  (*r)->state = memcpy((*r)->state, s->state, s->type->size);

  return GMCMC_SUCCESS;
}

/*!
 * Copies a random number generator into an existing random number generator.
 * If the source and destination RNGs have types with different sizes of state
 * this function will resize the destination RNG state vector.
 * 
 * @param [out] r  the destination RNG
 * @param [in]  s  the source RNG
 *
 * @return GMCMC_SUCCESS if the RNG was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the RNG state vector.
 */
gmcmc_error gmcmc_rng_copy(gmcmc_rng * r, const gmcmc_rng * s) {
  // If source == destination return now
  if (r == s)
    return GMCMC_SUCCESS;
  
  // If the RNGs have different sizes of state vectors
  if (r->type->size != s->type->size) {
    free(r->state);     // Free existing state vector

    // Allocate space for new state vector
    if ((r->state = malloc(s->type->size)) == NULL) {
      GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
      return GMCMC_ERROR_OUT_OF_MEMORY;
    }
  }

  // Copy the type and state vector
  r->type = memcpy(r->type, s->type, sizeof(gmcmc_rng_type));
  r->state = memcpy(r->state, s->state, s->type->size);

  return GMCMC_SUCCESS;
}

/*!
 * Destroys an RNG.
 *
 * @param [in] r  the RNG to destroy
 */
void gmcmc_rng_destroy(gmcmc_rng * r) {
  free(r->type);
  free(r->state);
  free(r);
}

/*!
 * (Re)seeds a random number generator.
 *
 * @param [in] r     the RNG to (re)seed
 * @param [in] seed  the new seed
 */
void gmcmc_rng_set(gmcmc_rng * r, uint64_t seed) {
  r->type->set(r->state, seed);
}

/*!
 * Generates a random integer uniformly distributed over the range [min, max].
 *
 * @param [in] r  the random number generator
 * @return a random integer.
 */
uint64_t gmcmc_rng_get(const gmcmc_rng * r) {
  return r->type->get(r->state);
}

/*!
 * Generates a random double precision real floating point value uniformly
 * distributed over the range (0, 1).
 *
 * @param [in] r  the random number generator
 * @return a random real.
 */
double gmcmc_rng_get_double(const gmcmc_rng * r) {
  return r->type->get_double(r->state);
}
