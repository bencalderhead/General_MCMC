#include "rng.h"
#include <stdlib.h>
#include <errno.h>

int rng_create(rng * r, const rng_type * type, uint64_t seed) {
  // Allocate space for the state vector
  if ((r->state = malloc(type->size)) == NULL)
    return ENOMEM;

  r->type = type;       // Store a pointer to the type of RNG created
  rng_set(r, seed);     // Seed the RNG to initialise the state vector

  return 0;
}

void rng_destroy(rng * r) {
  free(r->state);       // Free the state vector
}
