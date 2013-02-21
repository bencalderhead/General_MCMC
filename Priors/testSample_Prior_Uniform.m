function test_suite = testSample_Prior_Uniform
  initTestSuite;

function setup
  rng(0);

function samples = generate(lower, upper)
  samples = zeros(1, 100000);
  for i=1:100000
    samples(i) = Sample_Prior_Uniform([lower, upper]);
  end

% Test statistics
function testStatistics1
  lower = 0.0;
  upper = 1.0;
  X = generate(lower, upper);
  assertElementsAlmostEqual(0.5 * (lower + upper), mean(X), 'absolute', 0.000169);
  assertElementsAlmostEqual(0.08333333 * (upper - lower) * (upper - lower), var(X), 'absolute', 0.00003);

function testStatistics2
  lower = -1.0;
  upper = 1.0;
  X = generate(lower, upper);
  assertElementsAlmostEqual(0.5 * (lower + upper), mean(X), 'absolute', 0.00034);
  assertElementsAlmostEqual(0.08333333 * (upper - lower) * (upper - lower), var(X), 'absolute', 0.0002);

% Test parameter checking
%  function testUpperLower
%    assertExceptionThrown(@() Sample_Prior_Uniform([ 1.0, 0.0 ]), 'Uniform:Domain_Error');

% Test NaN propagation
function testNaNLower
  assertEqual(NaN, Sample_Prior_Uniform([ NaN, 1.0 ]));

function testNaNUpper
  assertEqual(NaN, Sample_Prior_Uniform([ 1.0, NaN ]));
