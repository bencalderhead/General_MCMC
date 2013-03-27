#ifndef GMCMC_ERROR_H
#define GMCMC_ERROR_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Geometric MCMC error type.
 */
typedef enum {
  GMCMC_SUCCESS = 0,                        /*!< No error */
  GMCMC_ERROR_OUT_OF_MEMORY = 1,            /*!< Out of memory */
  GMCMC_ERROR_INVALID_ARGUMENT = 2,         /*!< Invalid argument to a function */
} gmcmc_error;

/*!
 * Geometric MCMC error handler function type.
 */
typedef void (*gmcmc_error_handler_t)(gmcmc_error  /*!< [in] error that occured */,
                                      const char * /*!< [in] the expression that caused the error (may be NULL) */,
                                      const char * /*!< [in] the name of the calling function */,
                                      const char * /*!< [in] the name of the source file where the error occurred */,
                                      unsigned int /*!< [in] the line of the source file */);

/*!
 * The Geometric MCMC error handler.
 *
 * To get more detailed information about library errors when they occur set
 * this to an error handler function that will get called when an error occurs.
 *
 * To ignore errors set this to NULL.  Errors that occur will be returned from
 * the calling function and passed up the stack.
 *
 * The default error handler prints a message describing the error to the
 * console then continues.
 *
 * @see gmcmc_error_handler_t
 */
extern gmcmc_error_handler_t gmcmc_error_handler;

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
const char * gmcmc_error_string(gmcmc_error);

/*!
 * Error handler macro.
 *
 * Checks that the error handler is not NULL then invokes it.
 *
 * @see gmcmc_error_handler
 */
#define GMCMC_ERROR_HANDLER(error) \
  do { \
    if (gmcmc_error_handler != NULL) \
      gmcmc_error_handler(error, NULL, __func__, __FILE__, __LINE__); \
  } while (false)

/*!
 * Error checking macro.
 *
 * Wraps a function call or expression that evaluates to a Geometric MCMC error
 * code and invokes the error handler (if it is not NULL and) if the error code
 * is not GMCMC_SUCCESS then returns the error code.  The function, file and
 * line arguments of the error handler are populated using preprocessor macros.
 *
 * @see gmcmc_error_handler
 */
#define GMCMC_ERROR_CHECK(call) \
  do { \
    gmcmc_error error = (call); \
    if (error != GMCMC_SUCCESS) { \
      if (gmcmc_error_handler != NULL) \
        gmcmc_error_handler(error, #call, __func__, __FILE__, __LINE__); \
      return error; \
    } \
  } while (false)

#ifdef __cplusplus
}
#endif

#endif
