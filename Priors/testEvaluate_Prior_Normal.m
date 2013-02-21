function test_suite = testEvaluate_Prior_Normal
  initTestSuite;

%
% Normal distribution:
%
% Parameters:
%  mu
%  sigma > 0.0
%
% Support:
%  -Inf <= x <= +Inf
%
% http://en.wikipedia.org/wiki/Normal_distribution
%

% Test standard functionality
function testStandard
  Order = [ 0, 1, 2 ];
  Expected = [ 0.39894228, 0.0, -1.0 ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_Normal([ 0.0, 1.0 ], 0.0, Order(i)));
  end

% Test boundaries of support
function testNegInfX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, Inf, -1.0 ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_Normal([ 0.0, 1.0 ], -Inf, Order(i)));
  end

function testInfX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, -Inf, -1.0 ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_Normal([ 0.0, 1.0 ], Inf, Order(i)));
  end

% Test parameter checking
%  function testZeroSigma
%    Order = [ 0, 1, 2 ];
%    for i = 1:3
%      assertExceptionThrown(@() Evaluate_Prior_Normal([ 0.0, 0.0 ], 0.0, Order), 'Normal:Domain_Error');
%    end
%
%  function testOrder
%    assertExceptionThrown(@() Evaluate_Prior_Normal([ 1.0, 1.0 ], 1.0, 3), 'Prior:Domain_Error');

% Test Inf propagation
function testInfMu
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, Inf, -1.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Normal([ Inf, 1.0 ], 0.0, Order(i)));
  end

function testInfSigma
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(0.0, Evaluate_Prior_Normal([ 0.0, Inf ], 0.0, Order(i)));
  end

% Test NaN propagation
function testNaNX
  Order = [ 0, 1, 2 ];
  Expected = [ NaN, NaN, -1.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Normal([ 0.0, 1.0 ], NaN, Order(i)));
  end

function testNaNMu
  Order = [ 0, 1, 2 ];
  Expected = [ NaN, NaN, -1.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Normal([ NaN, 1.0 ], 0.0, Order(i)));
  end

function testNaNSigma
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(NaN, Evaluate_Prior_Normal([ 0.0, NaN ], 0.0, Order(i)));
  end
