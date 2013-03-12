#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include "gmcmc/error.h"
#include "gmcmc/rng.h"

// Test functions
static void test_rng_create();
static void test_rng_get();
static void test_rng_get_double();

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
  if ((suite = CU_add_suite("rng", NULL, NULL)) == NULL) {
    fprintf(stderr, "failed to create test suite: %s\n", CU_get_error_msg());
    return error;
  }
  
  // Add the tests to the suite
  if (CU_ADD_TEST(suite, test_rng_create) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_rng_get) == NULL) {
    fprintf(stderr, "failed to add test: %s\n", CU_get_error_msg());
    return error;
  }
  if (CU_ADD_TEST(suite, test_rng_get_double) == NULL) {
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

// Dummy RNG type
static void set(void * r, uint64_t s) {
  // Use seed as state
  uint64_t * state = (uint64_t *)r;
  *state = s;
}
static uint64_t get(void * r)      {
  // return state as uint64_t
  return *(uint64_t *)r;
}
static double get_double(void * r) {
  // return state cast as a double
  return (double)(*(uint64_t *)r);
}

static const gmcmc_rng_type type = { "test",
                                     0,
                                     UINT64_MAX,
                                     sizeof(uint64_t),
                                     set,
                                     get,
                                     get_double
                                   };

static void test_rng_create() {
  // Test failures for invalid RNG types
  gmcmc_rng * rng;
   
  // NULL name
  gmcmc_rng_type null_name = type;
  null_name.name = NULL;
  liberror = GMCMC_SUCCESS;     // Reset the global error status
  CU_ASSERT(gmcmc_rng_create(&rng, &null_name, 0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);  // Check that the error handler was called with the correct error code
   
  // min >= max
  gmcmc_rng_type min_max = type;
  min_max.min = UINT64_MAX;
  min_max.max = UINT64_MAX;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_rng_create(&rng, &min_max, 0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL seed function
  gmcmc_rng_type null_set = type;
  null_set.set = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_rng_create(&rng, &null_set, 0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL integer generator function
  gmcmc_rng_type null_get = type;
  null_get.get = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_rng_create(&rng, &null_get, 0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
   
  // NULL real generator function
  gmcmc_rng_type null_get_double = type;
  null_get_double.get_double = NULL;
  liberror = GMCMC_SUCCESS;
  CU_ASSERT(gmcmc_rng_create(&rng, &null_get_double, 0) == GMCMC_ERROR_INVALID_ARGUMENT);
  CU_ASSERT(liberror == GMCMC_ERROR_INVALID_ARGUMENT);
}

static void test_rng_get() {
  // Test integer generation function (returns seed)
  gmcmc_rng * rng;
  
  CU_ASSERT(gmcmc_rng_create(&rng, &type, 0) == GMCMC_SUCCESS);
  CU_ASSERT_EQUAL(gmcmc_rng_get(rng), 0);

  // Use seed function to re-seed RNG
  gmcmc_rng_set(rng, 5);
  CU_ASSERT_EQUAL(gmcmc_rng_get(rng), 5);

  gmcmc_rng_destroy(rng);
}

static void test_rng_get_double() {
  // Test real generation function (returns seed)
  gmcmc_rng * rng;
  
  CU_ASSERT(gmcmc_rng_create(&rng, &type, 0) == GMCMC_SUCCESS);
  CU_ASSERT_DOUBLE_EQUAL(gmcmc_rng_get_double(rng), 0.0, 0.0);

  // Use seed function to re-seed RNG
  gmcmc_rng_set(rng, 5);
  CU_ASSERT_DOUBLE_EQUAL(gmcmc_rng_get_double(rng), 5.0, 0.0);

  gmcmc_rng_destroy(rng);
}
