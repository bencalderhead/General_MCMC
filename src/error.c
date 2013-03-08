#include "mcmc/error.h"

#include <stdio.h>

static void default_error_handler(mcmc_error error, const char * call,
                                  const char * function, const char * file,
                                  unsigned int line) {
  if (call == NULL)
    fprintf(stderr, "MCMC error occurred in %s (%s:%u):\n"
                    "\t%s\n", call, function, file, line, mcmc_error_string(error));
  else
    fprintf(stderr, "MCMC error occurred when executing\n"
                    "  \"%s\"\n"
                    "in %s (%s:%u):\n"
                    "\t%s\n", call, function, file, line, mcmc_error_string(error));
}

mcmc_error_handler_t mcmc_error_handler = &default_error_handler;

const char * mcmc_error_string(mcmc_error error) {
  switch (error) {
    case MCMC_SUCCESS:                          return "No error";
    case MCMC_ERROR_OUT_OF_MEMORY:              return "Out of memory";
    case MCMC_ERROR_INVALID_ARGUMENT:           return "Invalid argument to a function";
    case MCMC_ERROR_PROPOSAL_OUTWITH_PRIOR:     return "Proposed parameter values lie outside prior";
    case MCMC_ERROR_NUMERICAL_INSTABILITY:      return "Possible numerical instability";
    case MCMC_ERROR_LIKELIHOOD:                 return "Calculation of the likelihood returned zero";
    case MCMC_ERROR_EIGENVALUES:                return "Failed to compute eigenvalues";
    case MCMC_ERROR_INVERSE:                    return "Failed to invert matrix";
    default:                                    return "Unknown error code";
  }
}
