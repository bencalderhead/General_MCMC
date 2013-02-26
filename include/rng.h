#ifndef RNG_H
#define RNG_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Random number generator algorithm.
 */
typedef struct {
  const char * name;    // Name of the algorithm
  uint64_t min, max;    // Minimum and maximum integer values generated
  size_t size;          // Size of state
  void (*set)(void *, uint64_t);  // Seed function
  uint64_t (*get)(void *);        // Integer generation function
  double (*get_double)(void *);   // Real generation function on (0,1)
} rng_type;

/**
 * Random number generator.
 */
typedef struct {
  const rng_type * type;        // The type of PRNG algorithm used
  void * state;                 // RNG state
} rng;

/**
 * Creates a new random number generator
 *
 * @param r     the created RNG is returned through this pointer
 * @param type  the RNG algorithm to use
 * @param seed  the seed to initialise the RNG state with
 *
 * @return 0 on success or ENOMEM if there is not enough memory to allocate a
 *         new RNG.
 */
int rng_create(rng **, const rng_type *, uint64_t);

/**
 * Destroys an RNG.
 *
 * @param r  the RNG to destroy
 */
void rng_destroy(rng *);

/**
 * (Re)seeds a random number generator.
 *
 * @param r     the RNG to (re)seed
 * @param seed  the new seed
 */
static inline void rng_set(rng * r, uint64_t seed) {
  r->type->set(r->state, seed);
}

/**
 * Generates a random integer uniformly distributed over the range [min, max].
 *
 * @param r  the random number generator
 * @return a random integer.
 */
static inline uint64_t rng_get(const rng * r) {
  return r->type->get(r->state);
}

/**
 * Generates a random double precision real floating point value uniformly
 * distributed over the range (0, 1).
 *
 * @param r  the random number generator
 * @return a random real.
 */
static inline double rng_get_double(const rng * r) {
  return r->type->get_double(r->state);
}

// Random number generator algorithms

/** Mersenne Twister 64-bit */
extern const rng_type * mt19337_64;

#ifdef __cplusplus
}
#endif

#endif
