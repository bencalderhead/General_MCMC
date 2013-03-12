#include "gmcmc/error.h"

#include <stdio.h>

/*!
 * Default error handler.
 * 
 * Prints a textual description of the error to the console.
 * 
 * @param [in] error  the error code
 * @param [in] call   the function call/expression that caused the error (may be
 *                      NULL)
 * @param [in] func   the name of the function where the error occurred
 * @param [in] file   the name of the source file where the error occurred
 * @param [in] line   the line number the error occurred on
 */
static void default_error_handler(gmcmc_error error, const char * call,
                                  const char * function, const char * file,
                                  unsigned int line) {
  if (call == NULL)
    fprintf(stderr, "Geometric MCMC error occurred in %s (%s:%u):\n"
                    "\t%s\n", function, file, line, gmcmc_error_string(error));
  else
    fprintf(stderr, "Geometric MCMC error occurred when executing\n"
                    "  \"%s\"\n"
                    "in %s (%s:%u):\n"
                    "\t%s\n", call, function, file, line, gmcmc_error_string(error));
}

/*!
 * The Geometric MCMC error handler.
 */
gmcmc_error_handler_t gmcmc_error_handler = &default_error_handler;

/*!
 * Error translation function.
 * 
 * Generates a textual description of Geometric MCMC error codes.
 * 
 * @param [in] error  the error code
 * 
 * @return a textual description of the error code.
 * 
 * @see gmcmc_error
 */
const char * gmcmc_error_string(gmcmc_error error) {
  switch (error) {
    case GMCMC_SUCCESS:                         return "No error";
    case GMCMC_ERROR_OUT_OF_MEMORY:             return "Out of memory";
    case GMCMC_ERROR_INVALID_ARGUMENT:          return "Invalid argument to a function";
    default:                                    return "Unknown error code";
  }
}
