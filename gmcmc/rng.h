#ifndef GMCMC_RNG_H
#define GMCMC_RNG_H

#include <gmcmc/error.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Random number generator algorithm.
 */
typedef struct __gmcmc_rng_type_st {
  const char * name;                    /**< Name of the algorithm */
  uint64_t min, max;                    /**< Minimum and maximum integer values generated */
  size_t size;                          /**< Size of state vector */
  void (*set)(void *, uint64_t);        /**< Seed function */
  uint64_t (*get)(void *);              /**< Integer generation function */
  double (*get_double)(void *);         /**< Real generation function on (0,1) */
} gmcmc_rng_type;

/*!
 * Random number generator.
 */
typedef struct __gmcmc_rng_st gmcmc_rng;

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
gmcmc_error gmcmc_rng_create(gmcmc_rng **, const gmcmc_rng_type *, uint64_t);

/*!
 * Creates a random number generator which is a duplicate of an existing random
 * number generator.
 *
 * @param [out] r    the random number generator to create
 * @param [in]  src  the random number generator to copy
 *
 * @return GMCMC_SUCCESS if the prior was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the RNG state vector.
 */
gmcmc_error gmcmc_rng_create_copy(gmcmc_rng **, const gmcmc_rng *);

/*!
 * Copies a random number generator into an existing random number generator.
 * If the source and destination RNGs have types with different sizes of state
 * this function will resize the destination RNG state vector.
 * 
 * @param [out] r    the destination RNG
 * @param [in]  src  the source RNG
 *
 * @return GMCMC_SUCCESS if the RNG was created successfully,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to allocate
 *         the RNG state vector.
 */
gmcmc_error gmcmc_rng_copy(gmcmc_rng *, const gmcmc_rng *);

/*!
 * Destroys an RNG.
 *
 * @param [in] r  the RNG to destroy
 */
void gmcmc_rng_destroy(gmcmc_rng *);

/*!
 * (Re)seeds a random number generator.
 *
 * @param [in] r     the RNG to (re)seed
 * @param [in] seed  the new seed
 */
void gmcmc_rng_set(gmcmc_rng *, uint64_t);

/*!
 * Generates a random integer uniformly distributed over the range [min, max].
 *
 * @param [in] r  the random number generator
 * @return a random integer.
 */
uint64_t gmcmc_rng_get(const gmcmc_rng *);

/*!
 * Generates a random double precision real floating point value uniformly
 * distributed over the range (0, 1).
 *
 * @param [in] r  the random number generator
 * @return a random real.
 */
double gmcmc_rng_get_double(const gmcmc_rng *);

// Random number generator algorithms

/*! Mersenne Twister 64-bit */
extern const gmcmc_rng_type * gmcmc_rng_mt19937_64;

/*! double-precision SIMD-oriented fast Mersenne Twister */
extern const gmcmc_rng_type * gmcmc_rng_dsfmt19937;

#ifdef __cplusplus
}
#endif

#endif
