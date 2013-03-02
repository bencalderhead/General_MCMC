#include "priors.h"

#include <math.h>

// sqrt(2.0 * M_PI)
#define M_SQRT2PI 2.50662827463100050241

static bool validate(const double * params) {
  return params[1] > 0.0;       // Check that shape > 0.0
}

static double sample(const rng * restrict r, const double * restrict params) {
  // Numerical Recipes in C++, 3rd edition (section 7.3, page 369)

  double u, v, q;
  do {
    u = rng_get_double(r);
    v = 1.7156 * (rng_get_double(r) - 0.5);
    double x = u - 0.449871;
    double y = fabs(v) + 0.386595;
    q = x * x + y * (0.19600 * y - 0.25472 * x);
  } while (q > 0.27597 && (q > 0.27846 || v * v > -4.0 * log(u) * u * u));

  return exp(params[0] + params[1] * v / u);
}

static double evaluate(double x, const double * params) {
  double logscale = params[0], shape = params[1];
  return (1.0 / (x * shape * M_SQRT2PI)) *
  exp(-((log(x) - logscale) * (log(x) - logscale)) / (2.0 * shape * shape));
}

static double evaluate_1st_order(double x, const double * params) {
  double logscale = params[0], shape = params[1];
  return -(1.0 / x) - ((log(x) - logscale) / (shape * shape)) * (1.0 / x);
}

static double evaluate_2nd_order(double x, const double * params) {
  (void)x;
  double logscale = params[0], shape = params[1];
  return (1.0 / (x * x)) + ((log(x) - logscale) / (shape * shape)) *
         (1.0 / (x * x)) - (1.0 / (shape * shape)) * (1.0 / (x * x));
}

static const prior_type type = { "Lognormal", validate, sample, evaluate,
                                 evaluate_1st_order, evaluate_2nd_order, 2 };

const prior_type * lognormal_prior = &type;
