#include <gmcmc/priors.h>
#include <math.h>

// sqrt(2.0 * M_PI)
#define M_SQRT2PI 2.50662827463100050241

typedef struct {
  double logscale, shape;
} lognormal;

static bool init(void * params, va_list list) {
  lognormal * l = (lognormal *)params;
  l->logscale = va_arg(list, double);
  l->shape = va_arg(list, double);
  return (l->shape > 0.0);
}

static double sample(const gmcmc_rng * r, const void * params) {
  lognormal * l = (lognormal *)params;
  // Numerical Recipes in C++, 3rd edition (section 7.3, page 369)

  double u, v, q;
  do {
    u = gmcmc_rng_get_double(r);
    v = 1.7156 * (gmcmc_rng_get_double(r) - 0.5);
    double x = u - 0.449871;
    double y = fabs(v) + 0.386595;
    q = x * x + y * (0.19600 * y - 0.25472 * x);
  } while (q > 0.27597 && (q > 0.27846 || v * v > -4.0 * log(u) * u * u));

  return exp(l->logscale + l->shape * v / u);
}

static double evaluate(double x, const void * params) {
  lognormal * l = (lognormal *)params;
  return (1.0 / (x * l->shape * M_SQRT2PI)) *
  exp(-((log(x) - l->logscale) * (log(x) - l->logscale)) / (2.0 * l->shape * l->shape));
}

static double evaluate_1st_order(double x, const void * params) {
  lognormal * l = (lognormal *)params;
  return -(1.0 / x) - ((log(x) - l->logscale) / (l->shape * l->shape)) * (1.0 / x);
}

static double evaluate_2nd_order(double x, const void * params) {
  (void)x;
  lognormal * l = (lognormal *)params;
  return (1.0 / (x * x)) + ((log(x) - l->logscale) / (l->shape * l->shape)) *
         (1.0 / (x * x)) - (1.0 / (l->shape * l->shape)) * (1.0 / (x * x));
}

static const struct __gmcmc_prior_type_st type = { "Lognormal", init, sample, evaluate,
                                                   evaluate_1st_order, evaluate_2nd_order, sizeof(lognormal) };

const gmcmc_prior_type * gmcmc_prior_lognormal = &type;
