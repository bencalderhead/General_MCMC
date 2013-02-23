#ifndef LOGNORMAL_H
#define LOGNORMAL_H

#include "../rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * LogNormal distribution.
 *
 * Parameters:
 *   shape > 0
 *   log-scale
 *
 * Ref: http://en.wikipedia.org/wiki/Lognormal_distribution
 */
typedef struct __lognormal_distribution_st * lognormal_distribution;

int lognormal_create(lognormal_distribution *, double, double);
void lognormal_destroy(lognormal_distribution);

double lognormal_get_shape(const lognormal_distribution);
double lognormal_get_scale(const lognormal_distribution);

int lognormal_set_shape(lognormal_distribution, double);
int lognormal_set_scale(lognormal_distribution, double);

double lognormal_rand(const rng, const lognormal_distribution);

double lognormal_pdf(double, const lognormal_distribution);
double lognormal_1st_order_pdf(double, const lognormal_distribution);
double lognormal_2nd_order_pdf(double, const lognormal_distribution);

#ifdef __cplusplus
}
#endif

#endif
