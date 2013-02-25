#include "rng.h"
#include <stdint.h>

/*
 * Mersenne Twister for 64-bit machines.
 * 
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
 */

#define NN 312
#define MM 156
#define MATRIX_A 0xB5026F5AA96619E9ULL
#define UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define LM 0x7FFFFFFFULL /* Least significant 31 bits */

struct mt_state {
  uint64_t state[NN];
  int i;
};

/* initializes mt[NN] with a seed */
static void set(void * s, uint64_t seed) {
  struct mt_state * mt = (struct mt_state *)s;

  mt->state[0] = seed;
  for (mt->i = 1; mt->i < NN; mt->i++) 
    mt->state[mt->i] = (UINT64_C(6364136223846793005) * 
                        (mt->state[mt->i - 1] ^ (mt->state[mt->i - 1] >> 62))
                        + (uint64_t)mt->i);
}

/* generates a random number on [0, 2^64-1]-interval */
static uint64_t get(void * s) {
  static const uint64_t mag01[2]={UINT64_C(0), MATRIX_A};

  struct mt_state * mt = (struct mt_state *)s;

  if (mt->i >= NN) { /* generate NN words at one time */
    int i;
    for (i = 0; i < NN - MM; i++) {
      uint64_t x = (mt->state[i] & UM) | (mt->state[i + 1] & LM);
      mt->state[i] = mt->state[i + MM] ^ (x >> 1) ^ mag01[(int)(x & UINT64_C(1))];
    }
    for (;i < NN - 1; i++) {
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

/* generates a random number on [0,1)-real-interval */
static double get_double(void * s) {
  return (double)(get(s) >> 11) * (1.0 / 9007199254740992.0);
}

static const struct __rng_type_st mt19937_64_type = {
  "Mersenne Twister for 64-bit machines",
  UINT64_MAX,
  0,
  NN * sizeof(uint64_t),
  set,
  get,
  get_double,
};

const rng_type mt_19937_64 = (const rng_type)&mt19937_64_type;
