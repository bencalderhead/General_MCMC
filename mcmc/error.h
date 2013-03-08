#ifndef MCMC_ERROR_H
#define MCMC_ERROR_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  // General errors
  MCMC_SUCCESS = 0,                         // No error
  MCMC_ERROR_OUT_OF_MEMORY = 1,             // Out of memory
  MCMC_ERROR_INVALID_ARGUMENT = 2,          // Invalid argument to a function

  MCMC_ERROR_PROPOSAL_OUTWITH_PRIOR = 4,    // Proposed parameter values lie outside prior
  MCMC_ERROR_NUMERICAL_INSTABILITY = 5,     // Possible numerical instability
  MCMC_ERROR_LIKELIHOOD = 6,                // Calculation of the likelihood returned zero
  MCMC_ERROR_EIGENVALUES = 7,               // Failed to compute eigenvalues
  MCMC_ERROR_INVERSE = 8                    // Failed to invert matrix
} mcmc_error;

typedef void (*mcmc_error_handler_t)(mcmc_error, const char *, const char *, const char *, unsigned int);

extern mcmc_error_handler_t mcmc_error_handler;

const char * mcmc_error_string(mcmc_error);

#define MCMC_ERROR_HANDLER(error) \
  do { \
    if (mcmc_error_handler != NULL) \
      mcmc_error_handler(error, NULL, __func__, __FILE__, __LINE__); \
  } while (false)

#define MCMC_ERROR_CHECK(call) \
  do { \
    mcmc_error error = (call); \
    if (error != MCMC_SUCCESS) { \
      if (mcmc_error_handler != NULL) \
        mcmc_error_handler(error, #call, __func__, __FILE__, __LINE__); \
      return error; \
    } \
  } while (false)

#ifdef __cplusplus
}
#endif

#endif
