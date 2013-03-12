#include "gmcmc/priors.h"

#include <math.h>

// sqrt(2.0 * M_PI)
#define M_SQRT2PI 2.50662827463100050241

typedef struct {
  double mean, stddev;
} normal;

static bool init(void * params, va_list list) {
  normal * n = (normal *)params;
  n->mean = va_arg(list, double);
  n->stddev = va_arg(list, double);
  return (n->stddev > 0.0);
}

static double sample(const gmcmc_rng * r, const void * params) {
  normal * n = (normal *)params;
  // Numerical Recipes in C++, 3rd edition (section 7.3, page 369)

  double u, v, q;
  do {
    u = gmcmc_rng_get_double(r);
    v = 1.7156 * (gmcmc_rng_get_double(r) - 0.5);
    double x = u - 0.449871;
    double y = fabs(v) + 0.386595;
    q = x * x + y * (0.19600 * y - 0.25472 * x);
  } while (q > 0.27597 && (q > 0.27846 || v * v > -4.0 * log(u) * u * u));

  return n->mean + n->stddev * v / u;
}

static double evaluate(double x, const void * params) {
  normal * n = (normal *)params;
  return 1.0 / (n->stddev * M_SQRT2PI) *
  exp(-((x - n->mean) * (x - n->mean)) / (2.0 * n->stddev * n->stddev));
}

static double evaluate_1st_order(double x, const void * params) {
  normal * n = (normal *)params;
  return -(x - n->mean) / (n->stddev * n->stddev);
}

static double evaluate_2nd_order(double x, const void * params) {
  (void)x;
  normal * n = (normal *)params;
  return -1.0 / (n->stddev * n->stddev);
}

static const struct __gmcmc_prior_type_st type = { "Normal", init, sample, evaluate,
                                                   evaluate_1st_order, evaluate_2nd_order, sizeof(normal) };

const gmcmc_prior_type * gmcmc_prior_normal = &type;
