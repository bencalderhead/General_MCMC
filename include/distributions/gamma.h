#ifndef GAMMA_H
#define GAMMA_H

#include <math.h>
#include "rng.h"
#include "normal.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Gamma distribution.
 *
 * Parameters:
 *   alpha > 0
 *   beta > 0
 *
 * http://en.wikipedia.org/wiki/Gamma_distribution
 */
typedef struct {
  double alpha, beta;
} gamma;

static double gammln(double xx) {
  // Numerical Recipes in C++, 2nd edition (section 6.1, page 219)
  static const double cof[6] = { 76.18009172947146,      -86.50532032941677,
                                 24.01409824083091,       -1.231739572450155,
                                  0.1208650973866179e-2,  -0.5395239384953e-5 };

  double y = xx, x = xx;
  double tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  double ser = 1.000000000190015;
  for (int j = 0; j < 6; j++)
    ser += cof[j] / ++y;
  return -tmp + log(M_SQRT2PI * ser / x);
}

/**
 * Generates a gamma distributed double precision floating point variable.
 *
 * @param r  a random number generator
 * @param g  distribution parameters
 *
 * @return a gamma distributed double precision floating point variable.
 */
static inline double gamma_rand(const rng * r, const gamma * g) {
  // Numerical Recipes in C++, 3rd edition (section 7.3, page 370)
  const normal n = { 0.0, 1.0 };
  const double alpha = (g->alpha < 1.0) ? g->alpha + 1.0 : g->alpha;
  const double a1 = alpha - 1.0 / 3.0;
  const double a2 = 1.0 / sqrt(9.0 * a1);

  double u, v, x;
  do {
    do {
      x = normal_rand(r, &n);
      v = 1.0 + a2 * x;
    } while (v <= 0.0);

    v = v * v * v;
    u = rng_get_double(r);
  } while (u > 1.0 - 0.331 * x * x * x * x &&
           log(u) > 0.5 * x * x + a1 * (1.0 - v + log(v)));

  if (alpha == g->alpha)
    return a1 * v / g->beta;
  else {
    do { u = rng_get_double(r); } while (u == 0.0);
    return pow(u, 1.0 / g->alpha) * a1 * v / g->beta;
  }
}

/**
 * Calculates the value of the probability density function for the specified
 * gamma distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param g  distribution parameters
 *
 * @return the value of the PDF at point x.
 */
static inline double gamma_pdf(double x, const gamma * g) {
  // Numerical Recipes in C++, 3rd edition (section 6.14, page 331)
  if (x <= 0.0)
    return 0.0;
  const double fac = g->alpha * log(g->beta) - gammln(g->alpha);
  return exp(-g->beta * x + (g->alpha - 1.0) * log(x) + fac);
}

/**
 * Calculates the value of the first derivative of the probability density
 * function for the specified gamma distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param g  distribution parameters
 *
 * @return the value of the first derivative of the PDF at point x.
 */
double gamma_1st_order_pdf(double x, const gamma * g) {
  return (x <= 0.0) ? -HUGE_VAL : (g->alpha - 1.0) / x - 1.0 / g->beta;
}

/**
 * Calculates the value of the second derivative of the probability density
 * function for the specified gamma distribution at point x.
 *
 * @param x  the point to evaluate the PDF at
 * @param g  distribution parameters
 *
 * @return the value of the second derivative of the PDF at point x.
 */
double gamma_2nd_order_pdf(double x, const gamma * g) {
  return (x <= 0.0) ? -HUGE_VAL : -(g->alpha - 1.0) / (x * x);
}

#ifdef __cplusplus
}
#endif

#endif
