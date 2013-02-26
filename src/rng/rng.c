#include "rng.h"
#include <stdlib.h>
#include <errno.h>

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
int rng_create(rng ** r, const rng_type * type, uint64_t seed) {
  if ((*r = malloc(sizeof(rng))) == NULL)
    return ENOMEM;

  if (((*r)->state = malloc(type->size)) == NULL) {
    free(*r);
    return ENOMEM;
  }

  (*r)->type = type;
  type->set((*r)->state, seed);

  return 0;
}

/**
 * Destroys an RNG.
 *
 * @param r  the RNG to destroy
 */
void rng_destroy(rng * r) {
  free(r->state);
  free(r);
}
