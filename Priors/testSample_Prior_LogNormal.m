function test_suite = testSample_Prior_LogNormal
  initTestSuite;

function setup
  rng(0);

function samples = generate(mu, sigma)
  samples = zeros(1, 100000);
  for i=1:100000
    samples(i) = Sample_Prior_LogNormal([mu, sigma]);
  end

% Test statistics
function testStatistics1
  mu = 1.0;
  sigma = 1.0;
  X = generate(mu, sigma);
  assertElementsAlmostEqual(exp(mu + sigma/2), mean(X), 'absolute', 0.02);
  assertElementsAlmostEqual((exp(sigma) - 1) * exp(2 * mu + sigma), var(X), 'absolute', 0.03);

function testStatistics2
  mu = 1.0;
  sigma = 0.01;
  X = generate(mu, sigma);
  assertElementsAlmostEqual(exp(mu + sigma/2), mean(X), 'absolute', 0.02);
  assertElementsAlmostEqual((exp(sigma) - 1) * exp(2 * mu + sigma), var(X), 'absolute', 0.08);

% Test parameter checking
%  function testInfMean
%    assertEqual(Inf, Sample_Prior_LogNormal([ Inf, 1.0 ]));
%
%  function testZeroSigma
%    assertEqual(1.0, Sample_Prior_LogNormal([ 0.0, 0.0 ]));
%
%  function testInfSigma
%    assertEqual(Inf, Sample_Prior_LogNormal([ 0.0, Inf ]));

% Test NaN propagation
function testNaNMean
  assertEqual(NaN, Sample_Prior_LogNormal([ NaN, 1.0 ]));

function testNaNSigma
  assertEqual(NaN, Sample_Prior_LogNormal([ 1.0, NaN ]));
