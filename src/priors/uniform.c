#include "gmcmc/priors.h"

#include <math.h>

typedef struct {
  double lower, upper;
} uniform;

static bool init(void * params, va_list list) {
  uniform * u = (uniform *)params;
  u->lower = va_arg(list, double);
  u->upper = va_arg(list, double);
  return (u->lower < u->upper);
}

static double sample(const gmcmc_rng * r, const void * params) {
  uniform * u = (uniform *)params;
  return u->lower + gmcmc_rng_get_double(r) * (u->upper - u->lower);
}

static double evaluate(double x, const void * params) {
  uniform * u = (uniform *)params;
  return (x <= u->lower || x >= u->upper) ? 0.0 : 1.0 / (u->upper - u->lower);
}

static double evaluate_1st_order(double x, const void * params) {
  uniform * u = (uniform *)params;
  return (x <= u->lower || x >= u->upper) ? -HUGE_VAL : 0.0;
}

static double evaluate_2nd_order(double x, const void * params) {
  uniform * u = (uniform *)params;
  return (x <= u->lower || x >= u->upper) ? -HUGE_VAL : 0.0;
}

static const struct __gmcmc_prior_type_st type = { "Uniform", init, sample, evaluate,
                                                   evaluate_1st_order, evaluate_2nd_order, sizeof(uniform) };

const gmcmc_prior_type * gmcmc_prior_uniform = &type;
