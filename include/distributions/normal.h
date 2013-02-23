#ifndef NORMAL_H
#define NORMAL_H

#include "../rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Normal distribution.
 *
 * Parameters:
 *   mean
 *   variance > 0
 *
 * Ref: http://en.wikipedia.org/wiki/Normal_distribution
 */
typedef struct __normal_distribution_st * normal_distribution;

int normal_create(normal_distribution *, double, double);
void normal_destroy(normal_distribution);

double normal_get_shape(const normal_distribution);
double normal_get_variance(const normal_distribution);

int normal_set_shape(normal_distribution, double);
int normal_set_variance(normal_distribution, double);

double normal_rand(const rng, const normal_distribution);

double normal_log_pdf(double, const normal_distribution);

double normal_pdf(double, const normal_distribution);
double normal_1st_order_pdf(double, const normal_distribution);
double normal_2nd_order_pdf(double, const normal_distribution);

#ifdef __cplusplus
}
#endif

#endif
