#include "rng.h"

#define NN 312
#define MM 156
#define MATRIX_A UINT64_C(0xB5026F5AA96619E9)
#define UM UINT64_C(0xFFFFFFFF80000000) /* Most significant 33 bits */
#define LM UINT64_C(0x7FFFFFFF) /* Least significant 31 bits */

/**
 * Mersenne Twister 64-bit.
 *
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
 */
typedef struct {
  uint64_t state[NN];   /* The array for the state vector */
  int i;
} mt64;

/**
 * Re-seeds the Mersenne Twister.
 *
 * @param rng    the RNG to re-seed
 * @param seed  the new seed
 */
static void set(void * rng, uint64_t seed) {
  mt64 * mt = (mt64 *)rng;

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
 * @param rng  the RNG
 * @return the random integer.
 */
static uint64_t get(void * rng) {
  mt64 * mt = (mt64 *)rng;

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
 * (0, 1).
 *
 * @param mt  the PRNG
 * @return the random double precision floating point value.
 */
static double get_double(void * rng) {
  mt64 * mt = (mt64 *)rng;
  return ((double)get(mt) + 1.0) / ((double)UINT64_MAX + 2.0);
}

static const rng_type mt64_type = {
  "Mersenne Twister 64-bit", 0, UINT64_MAX, sizeof(mt64), set, get, get_double
};

const rng_type * mt19937_64 = &mt64_type;