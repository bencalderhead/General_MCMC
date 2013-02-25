#include "rng/mt64.h"
#include <stdlib.h>
#include <errno.h>

#define NN 312
#define MM 156
#define MATRIX_A UINT64_C(0xB5026F5AA96619E9)
#define UM UINT64_C(0xFFFFFFFF80000000) /* Most significant 33 bits */
#define LM UINT64_C(0x7FFFFFFF) /* Least significant 31 bits */

/**
 * Mersenne Twister 64-bit version.
 *
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
 */
struct __mt64_st {
  uint64_t state[NN];   /* The array for the state vector */
  int i;
};

/**
 * Creates a new Mersenne Twister PRNG and initialises it with a seed.
 *
 * @param mt    the created PRNG is returned through this pointer
 * @param seed  the seed to use to initialise the PRNG state
 * @return 0 on success, or ENOMEM if there is not enough memory available to
 * create another Mersenne Twister PRNG.
 */
int mt64_create(mt64 * mt, uint64_t seed) {
  if ((*mt = malloc(sizeof(struct __mt64_st))) == NULL)
    return ENOMEM;

  mt64_set(*mt, seed);

  return 0;
}

/**
 * Destroys a Mersenne Twister PRNG, freeing any memory used.
 *
 * @param mt  the PRNG to destroy
 */
void mt64_destroy(mt64 mt) {
  free(mt);
}

/**
 * Re-seeds the Mersenne Twister.
 *
 * @param mt    the PRNG to re-seed
 * @param seed  the new seed
 */
void mt64_set(const mt64 mt, uint64_t seed) {
  mt->state[0] = seed;
  for (int i = 1; i < NN; i++)
    mt->state[i] =  (UINT64_C(6364136223846793005) *
                     (mt->state[i - 1] ^ (mt->state[i - 1] >> 62))
                     + (uint64_t)i);
  mt->i = 0;
}

/**
 * Generates a 64-bit unsigned integer uniformly distributed over [0, 2^64-1].
 *
 * @param mt  the PRNG
 * @return the random integer.
 */
uint64_t mt64_get(const mt64 mt) {
  static const uint64_t mag01[2] = { UINT64_C(0), MATRIX_A };

  if (mt->i >= NN) { /* generate NN words at one time */
    for (int i = 0; i < NN - MM; i++) {
      uint64_t x = (mt->state[i] & UM) | (mt->state[i + 1] & LM);
      mt->state[i] = mt->state[i + MM] ^ (x >> 1) ^ mag01[(int)(x & UINT64_C(1))];
    }

    for (int i = NN - MM; i < NN - 1; i++) {
      uint64_t x = (mt->state[i] & UM) | (mt->state[i + 1] & LM);
      mt->state[i] = mt->state[i + (MM - NN)] ^ (x >> 1) ^ mag01[(int)(x & UINT64_C(1))];
    }

    uint64_t x = (mt->state[NN - 1] & UM) | (mt->state[0] & LM);
    mt->state[NN - 1] = mt->state[MM - 1] ^ (x >> 1) ^ mag01[(int)(x & UINT64_C(1))];

    mt->i = 0;
  }

  uint64_t x = mt->state[mt->i++];

  x ^= (x >> 29) & UINT64_C(0x5555555555555555);
  x ^= (x << 17) & UINT64_C(0x71D67FFFEDA60000);
  x ^= (x << 37) & UINT64_C(0xFFF7EEE000000000);
  x ^= (x >> 43);

  return x;
}

/**
 * Generates a double precision floating point value uniformly distributed over
 * [0, 1).
 *
 * @param mt  the PRNG
 * @return the random double precision floating point value.
 */
double mt64_get_double(const mt64 mt) {
  return (double)mt64_get(mt) / ((double)UINT64_MAX + 1.0);
}
