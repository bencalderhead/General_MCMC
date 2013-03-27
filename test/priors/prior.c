#include <gmcmc/error.h>
#include <gmcmc/priors.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

// Test functions
static void test_prior_create();
static void test_prior_sample();
static void test_prior_evaluate();
static void test_prior_evaluate_1st_order();
static void test_prior_evaluate_2nd_order();

// Global library error status
static gmcmc_error liberror;

// Custom error handler
static void cu_error_handler(gmcmc_error error, const char * call,
                             const char * func, const char * file,
                             unsigned int line) {
  (void)call;
  (void)func;
  (void)file;
  (void)line;
  // Copy the library error into the global error status
  liberror = error;
}

int main() {
  // Set the library error handler to the custom one
  gmcmc_error_handler = cu_error_handler;

  // Initialise the CUnit test registry
  CU_ErrorCode error;
  if ((error = CU_initialize_registry()) != CUE_SUCCESS) {
    fprintf(stderr, "failed to initialise test registry: %s\n", CU_get_error_msg());
    return error;
  }
  
  // Create a test suite within the registry
  CU_pSuite suite;
  if ((suite = CU_add_suite("prior", NULL, NULL)) == NULL) {
    fprintf(stderr, "failed to create test suite: %s\n", CU_get_error_msg());
    return error;
  }
  
  // Add the tests to the suite
  if (CU_ADD_TEST(suite, test_prior_create) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_prior_sample) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_prior_evaluate) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_prior_evaluate_1st_order) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_prior_evaluate_2nd_order) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }

  // Run the test suites using the CUnit basic interface
  if ((error = CU_basic_run_tests()) != CUE_SUCCESS) {
    fprintf(stderr, "failed to run tests: %s\n", CU_get_error_msg());
    return error;
  }

  // Display any failures (plus hack for newline afterwards)
  CU_basic_show_failures(CU_get_failure_list());
  if (CU_get_number_of_failure_records() > 0) printf("\n");

  // Get the number of tests that failed
  unsigned int failures = CU_get_number_of_tests_failed();

  // Cleanup the test registry
  CU_cleanup_registry();

  // Return the number of test failures
  return (int)failures;
}

// Dummy prior type
static bool init(void * v, va_list ap) {
  double * d = (double *)v;
  d[0] = va_arg(ap, double);
  d[1] = va_arg(ap, double);
  d[2] = va_arg(ap, double);
  d[3] = va_arg(ap, double);
  return (d[0] < d[1] && d[1] < d[2] && d[2] < d[3]);
}
static double sample(const gmcmc_rng * r, const void * v) {
  (void)r;
  return ((double *)v)[0];
}
static double evaluate(double x, const void * v) { (void)x; return ((double *)v)[1]; }
static double evaluate_1st_order(double x, const void * v) { (void)x; return ((double *)v)[2]; }
static double evaluate_2nd_order(double x, const void * v) { (void)x; return ((double *)v)[3]; }

static const gmcmc_prior_type type = { "test",
                                       init,
                                       sample,
                                       evaluate,
                                       evaluate_1st_order,
                                       evaluate_2nd_order,
                                       4 * sizeof(double)
                                     };
                                       

static void test_prior_create() {
  // Test failures for invalid prior types
  gmcmc_prior * prior;
   
  // NULL name
  gmcmc_prior_type null_name = type;
  null_name.name = NULL;
  liberror = GMCMC_SUCCESS;     // Reset the global error status
  CU_ASSERT(gmcmc_prior_create(&prior, &null_name, 0.0, 1.0, 2.0, 3.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);  // Check that the error handler was called with the correct error code
   
  // NULL init function
  gmcmc_prior_type null_init = type;
  null_init.init = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_prior_create(&prior, &null_init, 0.0, 1.0, 2.0, 3.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL sample function
  gmcmc_prior_type null_sample = type;
  null_sample.sample = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_prior_create(&prior, &null_sample, 0.0, 1.0, 2.0, 3.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL evaluate function
  gmcmc_prior_type null_eval = type;
  null_eval.evaluate = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_prior_create(&prior, &null_eval, 0.0, 1.0, 2.0, 3.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL 1st order evaluate function
  gmcmc_prior_type null_eval1 = type;
  null_eval1.evaluate_1st_order = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_prior_create(&prior, &null_eval1, 0.0, 1.0, 2.0, 3.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL 2nd order evaluate function
  gmcmc_prior_type null_eval2 = type;
  null_eval2.evaluate_2nd_order = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_prior_create(&prior, &null_eval2, 0.0, 1.0, 2.0, 3.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL invalid arguments
  gmcmc_prior_type inv_args = type;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_prior_create(&prior, &inv_args, 3.0, 2.0, 1.0, 0.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
}

static void test_prior_sample() {
  // Test prior sample function (returns 1st parameter)
  gmcmc_prior * prior;
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 0.0, 1.0, 2.0, 3.0) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_sample(prior, NULL), 0.0);
  gmcmc_prior_destroy(prior);
  
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 1.0, 1.1, 1.2, 1.3) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_sample(prior, NULL), 1.0);
  gmcmc_prior_destroy(prior);
}

static void test_prior_evaluate() {
  // Test evaluate function (returns 2nd parameter)
  gmcmc_prior * prior;
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 0.0, 1.0, 2.0, 3.0) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_evaluate(prior, 0.0), 1.0);
  gmcmc_prior_destroy(prior);
  
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 1.0, 1.1, 1.2, 1.3) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_evaluate(prior, 0.0), 1.1);
  gmcmc_prior_destroy(prior);
}

static void test_prior_evaluate_1st_order() {
  // Test evaluate function (returns 3rd parameter)
  gmcmc_prior * prior;
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 0.0, 1.0, 2.0, 3.0) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_evaluate_1st_order(prior, 0.0), 2.0);
  gmcmc_prior_destroy(prior);
  
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 1.0, 1.1, 1.2, 1.3) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_evaluate_1st_order(prior, 0.0), 1.2);
  gmcmc_prior_destroy(prior);
}

static void test_prior_evaluate_2nd_order() {
  // Test evaluate function (returns 4th parameter)
  gmcmc_prior * prior;
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 0.0, 1.0, 2.0, 3.0) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_evaluate_2nd_order(prior, 0.0), 3.0);
  gmcmc_prior_destroy(prior);
  
  CU_ASSERT(gmcmc_prior_create(&prior, &type, 1.0, 1.1, 1.2, 1.3) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_prior_evaluate_2nd_order(prior, 0.0), 1.3);
  gmcmc_prior_destroy(prior);
}
