function test_suite = testSample_Prior_Gamma
  initTestSuite;

function setup
  rng(0);

function samples = generate(shape, scale)
  samples = zeros(1, 100000);
  for i=1:100000
    samples(i) = Sample_Prior_Gamma([shape, scale]);
  end

% Test statistics
function testStatistics1
  shape = 1.0;
  scale = 1.0;
  X = generate(shape, scale);
  assertElementsAlmostEqual(shape * scale, mean(X), 'absolute', 0.001);
  assertElementsAlmostEqual(shape * scale * scale, var(X), 'absolute', 0.01);

function testStatistics2
  shape = 1.0;
  scale = 0.01;
  X = generate(shape, scale);
  assertElementsAlmostEqual(shape * scale, mean(X), 'absolute', 0.00001);
  assertElementsAlmostEqual(shape * scale * scale, var(X), 'absolute', 0.00001);

% Test parameter checking
%  function testZeroShape
%    assertExceptionThrown(@() Sample_Prior_Gamma([ 0.0, 1.0 ]), 'Gamma:Domain_Error');
%
%  function testInfShape
%    assertExceptionThrown(@() Sample_Prior_Gamma([ Inf, 1.0 ]), 'Gamma:Domain_Error');
%
%  function testZeroScale
%    assertExceptionThrown(@() Sample_Prior_Gamma([ 1.0, 0.0 ]), 'Gamma:Domain_Error');
%
%  function testInfScale
%    assertExceptionThrown(@() Sample_Prior_Gamma([ 1.0, Inf ]), 'Gamma:Domain_Error');

% Test NaN propagation
function testNaNShape
  assertEqual(NaN, Sample_Prior_Gamma([ NaN, 1.0 ]));

function testNaNScale
  assertEqual(NaN, Sample_Prior_Gamma([ 1.0, NaN ]));
