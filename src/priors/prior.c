#include "priors.h"

#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>

int prior_create(prior * restrict p, const prior_type * restrict type, ...) {
  // Allocate space for the parameter vector
  if ((p->params = malloc((size_t)type->n * sizeof(double))) == NULL)
    return ENOMEM;

  p->type = type;       // Store a pointer to the type of prior created

  // Initialise the parameter vector
  va_list ap;
  va_start(ap, type);
  for (int i = 0; i < type->n; i++)
    p->params[i] = va_arg(ap, double);
  va_end(ap);

  // Check parameter values are valid
  if (!type->validate(p->params)) {
    free(p->params);
    return EINVAL;
  }

  return 0;
}

void prior_destroy(prior * p) {
  free(p->params);      // Free the parameter vector
}
