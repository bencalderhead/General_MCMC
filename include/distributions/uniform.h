#ifndef UNIFORM_H
#define UNIFORM_H

#include "../rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Uniform distribution.
 *
 * Parameters:
 *   lower < upper
 *
 * Ref: http://en.wikipedia.org/wiki/Uniform_distribution
 */
typedef struct __uniform_distribution_st * uniform_distribution;

int uniform_create(uniform_distribution *, double, double);
void uniform_destroy(uniform_distribution);

double uniform_get_lower(const uniform_distribution);
double uniform_get_upper(const uniform_distribution);

int uniform_set_lower(uniform_distribution, double);
int uniform_set_upper(uniform_distribution, double);

double uniform_rand(const rng, const uniform_distribution);

double uniform_pdf(double, const uniform_distribution);
double uniform_1st_order_pdf(double, const uniform_distribution);
double uniform_2nd_order_pdf(double, const uniform_distribution);

#ifdef __cplusplus
}
#endif

#endif
