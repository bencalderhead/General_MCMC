#include <gmcmc/error.h>
#include <gmcmc/priors.h>
#include <math.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

// Test functions
static int test_init();
static int test_cleanup();
static void test_params();
static void test_bounds();
static void test_evaluation();
static void test_statistics();

int main() {
  // Initialise the CUnit test registry
  CU_ErrorCode error;
  if ((error = CU_initialize_registry()) != CUE_SUCCESS) {
    fprintf(stderr, "failed to initialise test registry: %s\n", CU_get_error_msg());
    return error;
  }

  // Create a test suite within the registry
  CU_pSuite suite;
  if ((suite = CU_add_suite("gamma", test_init, test_cleanup)) == NULL) {
    fprintf(stderr, "failed to create test suite: %s\n", CU_get_error_msg());
    return error;
  }

  // Add the tests to the suite
  if (CU_ADD_TEST(suite, test_params) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_bounds) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_evaluation) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_statistics) == NULL) {
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

// Test fixture
static gmcmc_rng * rng;

static int test_init() {
  return gmcmc_rng_create(&rng, gmcmc_rng_mt19937_64, 0);
}

static int test_cleanup() {
  gmcmc_rng_destroy(rng);
  return 0;
}

// Test parameter checking
static void test_params() {
  // Set the error handler to NULL
  gmcmc_error_handler_t handler = gmcmc_error_handler;
  gmcmc_error_handler = NULL;

  gmcmc_prior * gamma;
  CU_ASSERT(gmcmc_prior_create(&gamma, gmcmc_prior_gamma, 1.0, 0.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(gmcmc_prior_create(&gamma, gmcmc_prior_gamma, 0.0, 1.0) == GMCMC_ERROR_INVALID_ARGUMENT);

  // Restore the error handler
  gmcmc_error_handler = handler;
}

// Test boundaries of support
static void test_bounds() {
  // Create the distribution
  gmcmc_prior * gamma;
  CU_ASSERT(gmcmc_prior_create(&gamma, gmcmc_prior_gamma, 1.0, 1.0) == GMCMC_SUCCESS);

  CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate(gamma, -1.e-18), 0.0, 0.0);
  CU_ASSERT_EQUAL(isinf(gmcmc_prior_evaluate_1st_order(gamma, -1.e-18)), -1);
  CU_ASSERT_EQUAL(isinf(gmcmc_prior_evaluate_2nd_order(gamma, -1.e-18)), -1);

  gmcmc_prior_destroy(gamma);
}

// Test evaluation
static void test_evaluation() {
  // Create the distribution
  gmcmc_prior * gamma;
  CU_ASSERT(gmcmc_prior_create(&gamma, gmcmc_prior_gamma, 1.0, 1.0) == GMCMC_SUCCESS);

  CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate(gamma, 1.0), 0.36787944, 1.e-8);
  CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate_1st_order(gamma, 1.0), -1.0, 0.0);
  CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate_2nd_order(gamma, 1.0), 0.0, 0.0);

  gmcmc_prior_destroy(gamma);
}

// Test statistics
static void test_statistics() {
  // Reseed the RNG
  gmcmc_rng_set(rng, 0);

  // Create the distribution
  gmcmc_prior * gamma;
  CU_ASSERT(gmcmc_prior_create(&gamma, gmcmc_prior_gamma, 1.0, 1.0) == GMCMC_SUCCESS);

  double x[100000], mean = 0.0, var = 0.0;
  for (size_t i = 0; i < 100000; i++) {
    x[i] = gmcmc_prior_sample(gamma, rng);

    // Knuth - The Art of Computer Programming (vol 2. 1998 p.232)
    double mprev = mean;
    mean += (x[i] - mean) / ((double)i + 1.0);
    var += (x[i] - mean) * (x[i] - mprev);
  }
  var /= 100000.0;

  // Test mean and variance
  CU_ASSERT_DOUBLE_EQUAL(mean, 1.0, 0.001);
  CU_ASSERT_DOUBLE_EQUAL(var, 1.0, 0.01);

  gmcmc_prior_destroy(gamma);
}
