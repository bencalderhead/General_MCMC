#include <gmcmc/rng.h>
#include <math.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

// Test functions
static int test_init();
static int test_cleanup();
static void test_rng_get();
static void test_rng_get_double();

static uint64_t mean(const uint64_t *, size_t);

int main() {
  // Initialise the CUnit test registry
  CU_ErrorCode error;
  if ((error = CU_initialize_registry()) != CUE_SUCCESS) {
    fprintf(stderr, "failed to initialise test registry: %s\n", CU_get_error_msg());
    return error;
  }
  
  // Create a test suite within the registry
  CU_pSuite suite;
  if ((suite = CU_add_suite("mt64", test_init, test_cleanup)) == NULL) {
    fprintf(stderr, "failed to create test suite: %s\n", CU_get_error_msg());
    return error;
  }
  
  // Add the tests to the suite
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

// RNG under test
static gmcmc_rng * rng;

static int test_init() {
  return gmcmc_rng_create(&rng, gmcmc_rng_mt19937_64, 0);
}

static int test_cleanup() {
  gmcmc_rng_destroy(rng);
  return 0;
}

static void test_rng_get() {
  // Reseed the RNG
  gmcmc_rng_set(rng, 0);

  // Test integer generation
  uint64_t x[100000];
  for (size_t i = 0; i < 100000; i++) {
    x[i] = gmcmc_rng_get(rng);

    // Test min/max
    CU_ASSERT(x[i] >= gmcmc_rng_mt19937_64->min);
    CU_ASSERT(x[i] <= gmcmc_rng_mt19937_64->max);
  }

  // Test mean
  uint64_t actual = mean(x, 100000);
  uint64_t expected = UINT64_MAX >> 1;
  uint64_t diff = (actual > expected) ? actual - expected : expected - actual;
  CU_ASSERT(diff <= 1e18);
}

static void test_rng_get_double() {
  // Reseed the RNG
  gmcmc_rng_set(rng, 0);

  // Test real generation
  double mean = 0.0, var = 0.0;
  for (size_t i = 0; i < 100000; i++) {
    double x = gmcmc_rng_get_double(rng);

    // Test min/max
    CU_ASSERT(x < 1.0);
    CU_ASSERT(x > 0.0);

    // Knuth - The Art of Computer Programming (vol 2. 1998 p.232)
    double mprev = mean;
    mean += (x - mean) / ((double)i + 1.0);
    var += (x - mean) * (x - mprev);
  }
  var /= 1.e5;

  // Test mean and variance
  CU_ASSERT_DOUBLE_EQUAL(mean, 0.5, 0.00169);           // 1/2 (a + b)
  CU_ASSERT_DOUBLE_EQUAL(var, 0.08333333, 0.0003);      // 1/12 (b - a)^2
}

// Utility functions

/*!
 * Calculates the mean of two 64-bit unsigned integers avoiding overflow.
 * 
 * @param [in] a  a 64-bit unsigned integer
 * @param [in] b  a 64-bit unsigned integer
 * 
 * @return the mean of a and b.
 */
static uint64_t mean2(uint64_t a, uint64_t b) {
  return (a & b) + ((a ^ b) >> 1);
}

/*!
 * Recursively calculates the mean of an array of 64-bit unsigned integers.
 * 
 * @param [in] x  an array of 64-bit unsigned integers
 * @param [in] l  the leftmost index
 * @param [in] r  the rightmost index (inclusive)
 * 
 * @return the mean of the elements in x.
 */
static uint64_t meanx(const uint64_t * x, size_t l, size_t r) {
  // The mean of one element is the element itself
  if (l >= r)
    return x[r];

  size_t m = l + ((r - l) >> 1);        // Midpoint
  return mean2(meanx(x, l, m), meanx(x, m + 1, r));
}

/*!
 * Calculates the mean of an array of 64-bit unsigned integers.
 * 
 * @param [in] x  an array of 64-bit unsigned integers
 * @param [in] n  the length of the array
 * 
 * @return the mean of the elements in x.
 */
uint64_t mean(const uint64_t * x, size_t n) {
  if (n == 0)
    return UINT64_C(0);
  return meanx(x, 0, n - 1);
}
