#include <gmcmc/rng.h>

#define DSFMT_MEXP 19937
/*-----------------
  BASIC DEFINITIONS
  -----------------*/
/* Mersenne Exponent. The period of the sequence
 *  is a multiple of 2^DSFMT_MEXP-1.
 * #define DSFMT_MEXP 19937 */
/** DSFMT generator has an internal state array of 128-bit integers,
 * and N is its size. */
#define DSFMT_N ((DSFMT_MEXP - 128) / 104 + 1)
/** N32 is the size of internal state array when regarded as an array
 * of 32-bit integers.*/
#define DSFMT_N32 (DSFMT_N * 4)
/** N64 is the size of internal state array when regarded as an array
 * of 64-bit integers.*/
#define DSFMT_N64 (DSFMT_N * 2)

#define DSFMT_LOW_MASK  UINT64_C(0x000FFFFFFFFFFFFF)
#define DSFMT_HIGH_CONST UINT64_C(0x3FF0000000000000)
#define DSFMT_SR        12

/* for sse2 */
#if defined(HAVE_SSE2)
  #define SSE2_SHUFF 0x1b
#elif defined(HAVE_ALTIVEC)
  #if defined(__APPLE__)  /* For OSX */
    #define ALTI_SR (vector unsigned char)(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4)
    #define ALTI_SR_PERM \
        (vector unsigned char)(15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14)
    #define ALTI_SR_MSK \
        (vector unsigned int)(0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU)
    #define ALTI_PERM \
        (vector unsigned char)(12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3)
  #else
    #define ALTI_SR      {4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4}
    #define ALTI_SR_PERM {15,0,1,2,3,4,5,6,15,8,9,10,11,12,13,14}
    #define ALTI_SR_MSK  {0x000fffffU,0xffffffffU,0x000fffffU,0xffffffffU}
    #define ALTI_PERM    {12,13,14,15,8,9,10,11,4,5,6,7,0,1,2,3}
  #endif
#endif

#define DSFMT_POS1      117
#define DSFMT_SL1       19
#define DSFMT_MSK1      UINT64_C(0x000ffafffffffb3f)
#define DSFMT_MSK2      UINT64_C(0x000ffdfffc90fffd)
#define DSFMT_MSK32_1   0x000ffaffU
#define DSFMT_MSK32_2   0xfffffb3fU
#define DSFMT_MSK32_3   0x000ffdffU
#define DSFMT_MSK32_4   0xfc90fffdU
#define DSFMT_FIX1      UINT64_C(0x90014964b32f4329)
#define DSFMT_FIX2      UINT64_C(0x3b8d12ac548a7c7a)
#define DSFMT_PCV1      UINT64_C(0x3d84e1ac0dc82880)
#define DSFMT_PCV2      UINT64_C(0x0000000000000001)
#define DSFMT_IDSTR     "dSFMT2-19937:117-19:ffafffffffb3f-ffdfffc90fffd"

/* PARAMETERS FOR ALTIVEC */
#if defined(__APPLE__)  /* For OSX */
    #define ALTI_SL1    (vector unsigned char)(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)
    #define ALTI_SL1_PERM \
        (vector unsigned char)(2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1)
    #define ALTI_SL1_MSK \
        (vector unsigned int)(0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U)
    #define ALTI_MSK    (vector unsigned int)(DSFMT_MSK32_1, \
                        DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4)
#else   /* For OTHER OSs(Linux?) */
    #define ALTI_SL1    {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3}
    #define ALTI_SL1_PERM \
        {2,3,4,5,6,7,30,30,10,11,12,13,14,15,0,1}
    #define ALTI_SL1_MSK \
        {0xffffffffU,0xfff80000U,0xffffffffU,0xfff80000U}
    #define ALTI_MSK \
        {DSFMT_MSK32_1, DSFMT_MSK32_2, DSFMT_MSK32_3, DSFMT_MSK32_4}
#endif

/*------------------------------------------
  128-bit SIMD like data type for standard C
  ------------------------------------------*/
#if defined(HAVE_ALTIVEC)
#  if !defined(__APPLE__)
#    include <altivec.h>
#  endif
/** 128-bit data structure */
typedef union W128_T {
    vector unsigned int s;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
} w128_t;

static inline void do_recursion(w128_t *r, w128_t *a, w128_t * b,
                                w128_t * lung) {
  const vector unsigned char sl1 = ALTI_SL1;
  const vector unsigned char sl1_perm = ALTI_SL1_PERM;
  const vector unsigned int sl1_msk = ALTI_SL1_MSK;
  const vector unsigned char sr1 = ALTI_SR;
  const vector unsigned char sr1_perm = ALTI_SR_PERM;
  const vector unsigned int sr1_msk = ALTI_SR_MSK;
  const vector unsigned char perm = ALTI_PERM;
  const vector unsigned int msk1 = ALTI_MSK;

  vector unsigned int z = a->s;
  vector unsigned int w = lung->s;
  vector unsigned int x = vec_perm(w, (vector unsigned int)perm, perm);
  vector unsigned int y = vec_perm(z, (vector unsigned int)sl1_perm, sl1_perm);
  y = vec_sll(y, sl1);
  y = vec_and(y, sl1_msk);
  w = vec_xor(x, b->s);
  w = vec_xor(w, y);
  x = vec_perm(w, (vector unsigned int)sr1_perm, sr1_perm);
  x = vec_srl(x, sr1);
  x = vec_and(x, sr1_msk);
  y = vec_and(w, msk1);
  z = vec_xor(z, y);
  r->s = vec_xor(z, x);
  lung->s = w;
}
#elif defined(HAVE_SSE2)
#  include <emmintrin.h>

/** 128-bit data structure */
typedef union W128_T {
    __m128i si;
    __m128d sd;
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
} w128_t;

union X128I_T {
    uint64_t u[2];
    __m128i  i128;
};

union X128D_T {
    double d[2];
    __m128d d128;
};

/** mask data for sse2 */
static const union X128I_T sse2_param_mask = { { DSFMT_MSK1, DSFMT_MSK2 } };

/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param d a 128-bit part of the internal state array (I/O)
 */
static inlinec void do_recursion(w128_t *r, w128_t *a, w128_t *b, w128_t *u) {
  __m128i x = a->si;
  __m128i z = _mm_slli_epi64(x, DSFMT_SL1);
  __m128i y = _mm_shuffle_epi32(u->si, SSE2_SHUFF);
  z = _mm_xor_si128(z, b->si);
  y = _mm_xor_si128(y, z);

  __m128i v = _mm_srli_epi64(y, DSFMT_SR);
  __m128i w = _mm_and_si128(y, sse2_param_mask.i128);
  v = _mm_xor_si128(v, x);
  v = _mm_xor_si128(v, w);
  r->si = v;
  u->si = y;
}
#else  /* standard C */
/** 128-bit data structure */
typedef union W128_T {
    uint64_t u[2];
    uint32_t u32[4];
    double d[2];
} w128_t;

/**
 * This function represents the recursion formula.
 * @param r output 128-bit
 * @param a a 128-bit part of the internal state array
 * @param b a 128-bit part of the internal state array
 * @param lung a 128-bit part of the internal state array (I/O)
 */
static inline void do_recursion(w128_t *r, w128_t *a, w128_t * b,
                                w128_t * lung) {
  uint64_t t0 = a->u[0];
  uint64_t t1 = a->u[1];
  uint64_t L0 = lung->u[0];
  uint64_t L1 = lung->u[1];
  lung->u[0] = (t0 << DSFMT_SL1) ^ (L1 >> 32) ^ (L1 << 32) ^ b->u[0];
  lung->u[1] = (t1 << DSFMT_SL1) ^ (L0 >> 32) ^ (L0 << 32) ^ b->u[1];
  r->u[0] = (lung->u[0] >> DSFMT_SR) ^ (lung->u[0] & DSFMT_MSK1) ^ t0;
  r->u[1] = (lung->u[1] >> DSFMT_SR) ^ (lung->u[1] & DSFMT_MSK2) ^ t1;
}
#endif

/**
 * Double precision SIMD-oriented fast Mersenne Twister.
 *
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/#dSFMT
 */
typedef struct {
  w128_t state[DSFMT_N + 1];
  int i;
} dsfmt_t;

/**
 * This function simulate a 32-bit array index overlapped to 64-bit
 * array of LITTLE ENDIAN in BIG ENDIAN machine.
 */
static inline int idxof(int i) {
#if defined(DSFMT_BIG_ENDIAN)
  return i ^ 1;
#else
  return i;
#endif
}

static void dsfmt_gen_rand_all(dsfmt_t *dsfmt) {
  w128_t lung = dsfmt->state[DSFMT_N];

  do_recursion(&dsfmt->state[0], &dsfmt->state[0],
                &dsfmt->state[DSFMT_POS1], &lung);

  int i;
  for (i = 1; i < DSFMT_N - DSFMT_POS1; i++)
    do_recursion(&dsfmt->state[i], &dsfmt->state[i],
                   &dsfmt->state[i + DSFMT_POS1], &lung);

  for (; i < DSFMT_N; i++)
    do_recursion(&dsfmt->state[i], &dsfmt->state[i],
                   &dsfmt->state[i + DSFMT_POS1 - DSFMT_N], &lung);

  dsfmt->state[DSFMT_N] = lung;
}

static void set(void * rng, uint64_t seed) {
  dsfmt_t * dsfmt = (dsfmt_t *)rng;

  uint32_t * psfmt32 = &dsfmt->state[0].u32[0];
  psfmt32[idxof(0)] = (uint32_t)seed;
  for (int i = 1; i < (DSFMT_N + 1) * 4; i++)
    psfmt32[idxof(i)] = UINT32_C(1812433253) *
      (psfmt32[idxof(i - 1)] ^ (psfmt32[idxof(i - 1)] >> 30)) + (uint32_t)i;

  uint64_t * psfmt64 = &dsfmt->state[0].u[0];
  for (int i = 0; i < DSFMT_N * 2; i++)
    psfmt64[i] = (psfmt64[i] & DSFMT_LOW_MASK) | DSFMT_HIGH_CONST;

  const uint64_t pcv[2] = { DSFMT_PCV1, DSFMT_PCV2 };
  const uint64_t tmp[] = { (dsfmt->state[DSFMT_N].u[0] ^ DSFMT_FIX1),
                           (dsfmt->state[DSFMT_N].u[1] ^ DSFMT_FIX2) };

  uint64_t inner = tmp[0] & pcv[0];
  inner ^= tmp[1] & pcv[1];
  for (int i = 32; i > 0; i >>= 1)
    inner ^= inner >> i;
  inner &= 1;
  // check OK
  if (inner == 1)
    return;

  // check NG, and modification
#if (DSFMT_PCV2 & 1) == 1
  dsfmt->state[DSFMT_N].u[1] ^= 1;
#else
  for (int i = 1; i >= 0; i--) {
    uint64_t work = 1;
    for (int j = 0; j < 64; j++) {
      if ((work & pcv[i]) != 0) {
        dsfmt->state[DSFMT_N].u[i] ^= work;
        return;
      }
      work = work << 1;
    }
  }
#endif

  dsfmt->i = DSFMT_N64;
}

static uint64_t get(void * rng) {
  dsfmt_t * dsfmt = (dsfmt_t *)rng;
  uint64_t * psfmt64 = &dsfmt->state[0].u[0];

  if (dsfmt->i >= DSFMT_N64) {
      dsfmt_gen_rand_all(dsfmt);
      dsfmt->i = 0;
  }

  return psfmt64[dsfmt->i++];
}

static double get_double(void * rng) {
  dsfmt_t * dsfmt = (dsfmt_t *)rng;
  double * dsfmt64 = &dsfmt->state[0].d[0];
  union {
    double d;
    uint64_t u;
  } r;

  if (dsfmt->i >= DSFMT_N64) {
    dsfmt_gen_rand_all(dsfmt);
    dsfmt->i = 0;
  }

  r.d = dsfmt64[dsfmt->i++];
  r.u |= 1;
  return r.d - 1.0;
}

static const struct __gmcmc_rng_type_st dsfmt_type = {
  "Double precision SIMD-oriented fast Mersenne Twister", 0, UINT64_MAX,
  sizeof(dsfmt_t), set, get, get_double
};

const gmcmc_rng_type * gmcmc_rng_dsfmt19937 = &dsfmt_type;
