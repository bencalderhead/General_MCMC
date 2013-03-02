#ifndef PRIORS_H
#define PRIORS_H

#include <stdbool.h>
#include "rng.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Prior type.
 */
typedef struct {
  const char * name;    // Name of the distribution
  bool (*validate)(const double *);     // Function to check parameter values
  double (*sample)(const rng *, const double *);        // Sample function
  double (*evaluate)(double, const double *);           // Evaluate function
  double (*evaluate_1st_order)(double, const double *); // 1st order derivative evaluate function
  double (*evaluate_2nd_order)(double, const double *); // 2nd order derivative evaluate function
  int n;                // Number of parameters
} prior_type;

/**
 * Prior.
 */
typedef struct {
  const prior_type * type;      // Type of prior
  double * params;              // Parameter vector
} prior;

/**
 * Creates a prior distribution.
 *
 * @param p       the prior distribution to create
 * @param type    the type of prior distribution to create
 * @param params  distribution parameters
 *
 * @return 0 on success, non-zero if there is not enough memory to allocate
 * space for the parameter vector or the parameter values are invalid.
 */
int prior_create(prior *, const prior_type *, ...);

/**
 * Destroys a prior distribution.
 *
 * @param p the prior distribution
 */
void prior_destroy(prior *);

/**
 * Generates a sample from a prior distribution.
 *
 * @param p  a prior distribution
 * @param r  the random number generator to use
 *
 * @return a sample from the prior
 */
static inline double prior_sample(const prior * p, const rng * r) {
  return p->type->sample(r, p->params);
}

/**
 * Evaluates a prior distribution probability density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
static inline double prior_evaluate(const prior * p, double x) {
  return p->type->evaluate(x, p->params);
}

/**
 * Evaluates the first order derivative of a prior distribution probability
 * density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
static inline double prior_evaluate_1st_order(const prior * p, double x) {
  return p->type->evaluate_1st_order(x, p->params);
}

/**
 * Evaluates the second order derivative of a prior distribution probability
 * density function.
 *
 * @param p  a prior distribution
 * @param x  the point at which to evaluate the probability density function
 *
 * @return the value of the prior probability density function at x.
 */
static inline double prior_evaluate_2nd_order(const prior * p, double x) {
  return p->type->evaluate_2nd_order(x, p->params);
}

// Predefined priors

/** A uniform prior over (a, b) */
extern const prior_type * uniform_prior;

/** A gamma prior with parameters (a, b) */
extern const prior_type * gamma_prior;

/** A normal prior with parameters for mean and standard deviation */
extern const prior_type * normal_prior;

/** A lognormal prior with parameters for log-mean and standard deviation */
extern const prior_type * lognormal_prior;

#ifdef __cplusplus
}
#endif

#endif
