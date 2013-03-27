#include <gmcmc/priors.h>
#include <math.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

// Test functions
static int test_init();
static int test_cleanup();
static void test_params();
static void test_bounds();
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
  if ((suite = CU_add_suite("uniform", test_init, test_cleanup)) == NULL) {
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

  gmcmc_prior * uniform;
  CU_ASSERT(gmcmc_prior_create(&uniform, gmcmc_prior_uniform, 1.0, 0.0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(gmcmc_prior_create(&uniform, gmcmc_prior_uniform, 1.0, 1.0) == GMCMC_ERROR_INVALID_ARGUMENT);

  // Restore the error handler
  gmcmc_error_handler = handler;
}

// Test boundaries of support
static void test_bounds() {
  // Create the distribution
  gmcmc_prior * uniform;
  CU_ASSERT(gmcmc_prior_create(&uniform, gmcmc_prior_uniform, 0.0, 1.0) == GMCMC_SUCCESS);

  // Test upper bound
  CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate(uniform, 1.0), 0.0, 0.0);
  CU_ASSERT_EQUAL(isinf(gmcmc_prior_evaluate_1st_order(uniform, 1.0)), -1);
  CU_ASSERT_EQUAL(isinf(gmcmc_prior_evaluate_2nd_order(uniform, 1.0)), -1);

  // Test lower bound
  CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate(uniform, 0.0), 0.0, 0.0);
  CU_ASSERT_EQUAL(isinf(gmcmc_prior_evaluate_1st_order(uniform, 0.0)), -1);
  CU_ASSERT_EQUAL(isinf(gmcmc_prior_evaluate_2nd_order(uniform, 0.0)), -1);

  gmcmc_prior_destroy(uniform);
}

// Test statistics
static void test_statistics() {
  // Reseed the RNG
  gmcmc_rng_set(rng, 0);

  // Create the distribution
  gmcmc_prior * uniform;
  CU_ASSERT(gmcmc_prior_create(&uniform, gmcmc_prior_uniform, 0.0, 1.0) == GMCMC_SUCCESS);

  double x[100000], mean = 0.0, var = 0.0;
  for (size_t i = 0; i < 100000; i++) {
    x[i] = gmcmc_prior_sample(uniform, rng);

    // Test min/max
    CU_ASSERT(x[i] > 0.0);
    CU_ASSERT(x[i] < 1.0);

    // Test evaluation
    CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate(uniform, x[i]), 1.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate_1st_order(uniform, x[i]), 0.0, 0.0);
    CU_ASSERT_DOUBLE_EQUAL(gmcmc_prior_evaluate_2nd_order(uniform, x[i]), 0.0, 0.0);

    // Knuth - The Art of Computer Programming (vol 2. 1998 p.232)
    double mprev = mean;
    mean += (x[i] - mean) / ((double)i + 1.0);
    var += (x[i] - mean) * (x[i] - mprev);
  }
  var /= 100000.0;

  // Test mean and variance
  CU_ASSERT_DOUBLE_EQUAL(mean, 0.5, 0.00113969);        // 1/2  (a + b)
  CU_ASSERT_DOUBLE_EQUAL(var, 0.08333333, 0.00011269);  // 1/12 (b - a)^2

  gmcmc_prior_destroy(uniform);
}
