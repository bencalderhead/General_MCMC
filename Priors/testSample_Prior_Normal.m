function test_suite = testSample_Prior_Normal
  initTestSuite;

function setup
  rng(0);

function samples = generate(mu, sigma)
  samples = zeros(1, 100000);
  for i=1:100000
    samples(i) = Sample_Prior_Normal([mu, sigma]);
  end

% Test statistics
function testStatistics1
  mu = 1.0;
  sigma = 1.0;
  X = generate(mu, sigma);
  assertElementsAlmostEqual(mu, mean(X), 'absolute', 0.0008);
  assertElementsAlmostEqual(sigma, var(X), 'absolute', 0.007);

function testStatistics2
  mu = 1.0;
  sigma = 0.01;
  X = generate(mu, sigma);
  assertElementsAlmostEqual(mu, mean(X), 'absolute', 0.000008);
  assertElementsAlmostEqual(sigma, var(X), 'absolute', 0.01);

% Test parameter checking
%  function testInfMean
%    assertEqual(Inf, Sample_Prior_Normal([ Inf, 1.0 ]));
%
%  function testZeroSigma
%    assertEqual(0.0, Sample_Prior_Normal([ 0.0, 0.0 ]));
%
%  function testInfSigma
%    assertEqual(Inf, Sample_Prior_Normal([ 0.0, Inf ]));

% Test NaN propagation
function testNaNMean
  assertEqual(NaN, Sample_Prior_Normal([ NaN, 1.0 ]));

function testNaNSigma
  assertEqual(NaN, Sample_Prior_Normal([ 1.0, NaN ]));
