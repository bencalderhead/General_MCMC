#include <gmcmc/priors.h>
#include <math.h>

// sqrt(2.0 * M_PI)
#define M_SQRT2PI 2.50662827463100050241

typedef struct {
  double alpha, beta;
} gamma;

static double gammln(double x) {
  // Numerical Recipes in C++, 2nd edition (section 6.1, page 219)
  static const double cof[6] = { 76.18009172947146,      -86.50532032941677,
                                 24.01409824083091,       -1.231739572450155,
                                  0.1208650973866179e-2,  -0.5395239384953e-5 };

  double y = x;
  double tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  double ser = 1.000000000190015;
  for (int j = 0; j < 6; j++)
    ser += cof[j] / ++y;
  return -tmp + log(M_SQRT2PI * ser / x);
}

static bool init(void * params, va_list list) {
  gamma * g = (gamma *)params;
  g->alpha = va_arg(list, double);
  g->beta = va_arg(list, double);
  return (g->alpha > 0.0) && (g->beta > 0.0);
}

static double sample(const gmcmc_rng * r, const void * params) {
  gamma * g = (gamma *)params;

  // Numerical Recipes in C++, 3rd edition (section 7.3, page 370)
  const double alpha = (g->alpha < 1.0) ? g->alpha + 1.0 : g->alpha;
  const double a1 = alpha - 1.0 / 3.0;
  const double a2 = 1.0 / sqrt(9.0 * a1);

  double u, v, x, q;
  do {
    do {
      // Numerical Recipes in C++, 3rd edition (section 7.3, page 369)
      do {
        u = gmcmc_rng_get_double(r);
        v = 1.7156 * (gmcmc_rng_get_double(r) - 0.5);
        x = u - 0.449871;
        double y = fabs(v) + 0.386595;
        q = x * x + y * (0.19600 * y - 0.25472 * x);
      } while (q > 0.27597 && (q > 0.27846 || v * v > -4.0 * log(u) * u * u));

      x = v / u;
      v = 1.0 + a2 * x;
    } while (v <= 0.0);

    v = v * v * v;
    u = gmcmc_rng_get_double(r);
  } while (u > 1.0 - 0.331 * x * x * x * x &&
           log(u) > 0.5 * x * x + a1 * (1.0 - v + log(v)));

  return (alpha == g->alpha) ? a1 * v / g->beta :
         pow(gmcmc_rng_get_double(r), 1.0 / g->alpha) * a1 * v / g->beta;
}

static double evaluate(double x, const void * params) {
  gamma * g = (gamma *)params;
  // Numerical Recipes in C++, 3rd edition (section 6.14, page 331)
  const double fac = g->alpha * log(g->beta) - gammln(g->alpha);
  return (x <= 0.0) ? 0.0 : exp(-g->beta * x + (g->alpha - 1.0) * log(x) + fac);
}

static double evaluate_1st_order(double x, const void * params) {
  gamma * g = (gamma *)params;
  return (x <= 0.0) ? -INFINITY : (g->alpha - 1.0) / x - 1.0 / g->beta;
}

static double evaluate_2nd_order(double x, const void * params) {
  gamma * g = (gamma *)params;
  return (x <= 0.0) ? -INFINITY : -(g->alpha - 1.0) / (x * x);
}

static const struct __gmcmc_prior_type_st type = { "Gamma", init, sample, evaluate,
                                                   evaluate_1st_order, evaluate_2nd_order, sizeof(gamma) };

const gmcmc_prior_type * gmcmc_prior_gamma = &type;
