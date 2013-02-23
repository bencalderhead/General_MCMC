#ifndef GAMMA_H
#define GAMMA_H

#include "../rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Gamma distribution.
 *
 * Parameters:
 *   shape > 0
 *   scale > 0
 *
 * Ref: http://en.wikipedia.org/wiki/Gamma_distribution
 */
typedef struct __gamma_distribution_st * gamma_distribution;

int gamma_create(gamma_distribution *, double, double);
void gamma_destroy(gamma_distribution);

double gamma_get_shape(const gamma_distribution);
double gamma_get_scale(const gamma_distribution);

int gamma_set_shape(gamma_distribution, double);
int gamma_set_scale(gamma_distribution, double);

double gamma_rand(const rng, const gamma_distribution);

double gamma_pdf(double, const gamma_distribution);
double gamma_1st_order_pdf(double, const gamma_distribution);
double gamma_2nd_order_pdf(double, const gamma_distribution);

#ifdef __cplusplus
}
#endif

#endif
