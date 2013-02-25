#include "rng.h"
#include <stdlib.h>
#include <errno.h>

int rng_create(rng * r, const rng_type t, uint64_t seed) {
  if ((*r = malloc(sizeof(struct __rng_st))) == NULL)
    return ENOMEM;
  
  if (((*r)->state = malloc(t->size)) == NULL) {
    free(*r);
    return ENOMEM;
  }
  
  (*r)->type = t;
  
  t->set((*r)->state, seed);
  
  return 0;
}

void rng_destroy(rng r) {
  free(r->state);
  free(r);
}
