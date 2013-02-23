#ifndef RNG_H
#define RNG_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct __rng_type_st {
  const char * name;
  unsigned long max, min;
  size_t size;
  void (*set)(void *, unsigned long);
  unsigned long (*get)(void *);
  double (*get_double)(void *);
} * rng_type;

typedef struct __rng_st {
  const rng_type type;
  void * state;
} * rng;

int rng_create(rng *, const rng_type);
void rng_destroy(rng);

static inline void rng_set(rng r, unsigned long seed) {
  r->set(r->state, seed);
}

static inline int rng_get(const rng r) { return r->get(r->state); }
static inline double rng_get_double(const rng r) { return r->get_double(r->state); }

const rng_type mt19937;

#ifdef __cplusplus
}
#endif

#endif
