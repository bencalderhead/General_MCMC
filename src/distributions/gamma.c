#include "distributions/gamma.h"
#include "rng.h"
#include <math.h>
#include <errno.h>

struct __gamma_distribution_st {
  double shape, scale;
};

int gamma_create(gamma_distribution * g, double shape, double scale) {
  if (shape <= 0.0 || scale <= 0.0)
    return EINVAL;

  if ((*g = malloc(sizeof(struct __gamma_distribution_st))) == NULL)
    return ENOMEM;

  (*g)->shape = shape;
  (*g)->scale = scale;

  return 0;
}

void gamma_destroy(gamma_distribution g) {
  free(g);
}

double gamma_get_shape(const gamma_distribution g) { return g->shape; }
double gamma_get_scale(const gamma_distribution g) { return g->scale; }

int gamma_set_shape(gamma_distribution g, double shape) {
  if (shape <= 0.0)
    return EINVAL;
  g->shape = shape;
  return 0;
}

int gamma_set_scale(gamma_distribution g, double scale) {
  if (scale <= 0.0)
    return EINVAL;
  g->scale = scale;
  return 0;
}

static inline double gamma_large(const rng r, const double shape) {
  /* Works only if a > 1, and is most efficient if a is large

    This algorithm, reported in Knuth, is attributed to Ahrens.  A
    faster one, we are told, can be found in: J. H. Ahrens and
    U. Dieter, Computing 12 (1974) 223-246.  */

  double sqa, x, y, v;
  sqa = sqrt(2.0 * shape - 1.0);
  do
    {
      do
        {
          y = tan(M_PI * rng_get_double(r));
          x = sqa * y + shape - 1.0;
        }
      while (x <= 0.0);
      v = rng_get_double(r);
    }
  while (v > (1.0 + y * y) * exp((shape - 1.0) * log(x / (shape - 1.0)) - sqa * y));

  return x;
}

static inline double gamma_frac(const rng r, const double shape) {
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;

  p = M_E / (shape + M_E);
  do {
    u = rng_get_double(r);
    v = 1.0 - rng_get_double(r);

    if (u < p) {
      x = exp((1.0 / shape) * log(v));
      q = exp(-x);
    }
    else {
      x = 1.0 - log(v);
      q = exp((shape - 1.0) * log(x));
    }
  } while (rng_get_double(r) >= q);

  return x;
}

static double gsl_ran_gamma_int(const rng * r, const unsigned int a) {
  if (a < 12)
    {
      unsigned int i;
      double prod = 1;

      for (i = 0; i < a; i++)
        {
          prod *= gsl_rng_uniform_pos (r);
        }

      /* Note: for 12 iterations we are safe against underflow, since
         the smallest positive random number is O(2^-32). This means
         the smallest possible product is 2^(-12*32) = 10^-116 which
         is within the range of double precision. */

      return -log (prod);
    }
  else
    {
      return gamma_large (r, (double) a);
    }
}

double gamma_rand(const rng r, const gamma_distribution g) {
  double shape, scale;
  if (g == NULL) {
    shape = 1.0;
    scale = 1.0;
  }
  else {
    shape = g->shape;
    scale = g->scale;
  }

  if (shape < 1.0) {
    double u = 1.0 - rng_get_double(r);
    return gamma_rand(r, 1.0 + shape, scale) * pow(u, 1.0 / shape);
  }

  double x, v, u;
  double d = a - 1.0 / 3.0;
  double c = (1.0 / 3.0) / sqrt(d);

  while (true) {
    double x, v;
    do {
      x = gsl_ran_gaussian_ziggurat(r, 1.0);
      v = 1.0 + c * x;
    } while (v <= 0.0);

    v = v * v * v;
    double u = 1.0 - gamma_rand(r);

    if (u < 1.0 - 0.0331 * x * x * x * x)
      break;

    if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v)))
      break;
  }

  return b * d * v;
}

double gamma_pdf(double, const gamma_distribution);
double gamma_1st_order_pdf(double, const gamma_distribution);
double gamma_2nd_order_pdf(double, const gamma_distribution);
