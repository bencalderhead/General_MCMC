#ifndef RNG_H
#define RNG_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct __rng_type_st {
  const char * name;
  uint64_t max, min;
  size_t size;
  void (*set)(void *, uint64_t);
  uint64_t (*get)(void *);
  double (*get_double)(void *);
} * rng_type;

typedef struct __rng_st {
  const struct __rng_type_st * type;
  void * state;
} * rng;

int rng_create(rng *, const rng_type, uint64_t);
void rng_destroy(rng);

static inline void rng_set(rng r, uint64_t seed) {
  r->type->set(r->state, seed);
}

static inline uint64_t rng_get(const rng r) { return r->type->get(r->state); }
static inline double rng_get_double(const rng r) { return r->type->get_double(r->state); }

extern const rng_type mt19937_64;

#ifdef __cplusplus
}
#endif

#endif
