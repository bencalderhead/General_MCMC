#include "priors.h"

#include <math.h>

static bool validate(const double * params) {
  return (params[0] < params[1]);       // Check that lower < upper
}

static double sample(const rng * restrict r, const double * restrict params) {
  return params[0] + rng_get_double(r) * (params[1] - params[0]);
}

static double evaluate(double x, const double * params) {
  return (x <= params[0] || x >= params[1]) ? 0.0 : 1.0 / (params[1] - params[0]);
}

static double evaluate_1st_order(double x, const double * params) {
  return (x <= params[0] || x >= params[1]) ? -HUGE_VAL : 0.0;
}

static double evaluate_2nd_order(double x, const double * params) {
  return (x <= params[0] || x >= params[1]) ? -HUGE_VAL : 0.0;
}

static const prior_type type = { "Uniform", validate, sample, evaluate,
                                 evaluate_1st_order, evaluate_2nd_order, 2 };

const prior_type * uniform_prior = &type;
