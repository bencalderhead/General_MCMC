#ifndef MT64_H
#define MT64_H

#include <stdint.h>

/**
 * Mersenne Twister 64-bit version.
 *
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
 */
typedef struct __mt64_st * mt64;

/**
 * Creates a new Mersenne Twister PRNG and initialises it with a seed.
 *
 * @param mt    the created PRNG is returned through this pointer
 * @param seed  the seed to use to initialise the PRNG state
 * @return 0 on success, or ENOMEM if there is not enough memory available to
 * create another Mersenne Twister PRNG.
 */
int mt64_create(mt64 *, uint64_t);

/**
 * Destroys a Mersenne Twister PRNG, freeing any memory used.
 *
 * @param mt  the PRNG to destroy
 */
void mt64_destroy(mt64);

/**
 * Re-seeds the Mersenne Twister.
 *
 * @param mt    the PRNG to re-seed
 * @param seed  the new seed
 */
void mt64_set(const mt64, uint64_t);

/**
 * Generates a 64-bit unsigned integer uniformly distributed over [0, 2^64-1].
 *
 * @param mt  the PRNG
 * @return the random integer.
 */
uint64_t mt64_get(const mt64);

/**
 * Generates a double precision floating point value uniformly distributed over
 * [0, 1).
 *
 * @param mt  the PRNG
 * @return the random double precision floating point value.
 */
double mt64_get_double(const mt64);

#endif
