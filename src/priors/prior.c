#include "priors.h"

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

int prior_create(prior * restrict p, const prior_type * restrict type, ...) {
  // Allocate space for the parameter vector
  if ((p->params = malloc(type->n)) == NULL)
    return ENOMEM;

  p->type = type;       // Store a pointer to the type of prior created

  // Initialise the parameter vector
  va_list list;
  va_start(list, type);
  bool valid = type->init(p->params, list);
  va_end(list);

  // Check parameter values are valid
  if (!valid) {
    free(p->params);
    return EINVAL;
  }

  return 0;
}

void prior_destroy(prior * p) {
  free(p->params);      // Free the parameter vector
}

int prior_copy(prior * restrict dest, const prior * restrict src) {
  // Allocate space for the parameter vector
  if ((dest->params = malloc(src->type->n)) == NULL)
    return ENOMEM;

  // Copy the type
  dest->type = src->type;

  // Copy the parameters
  dest->params = memcpy(dest->params, src->params, dest->type->n);

  return 0;
}
