function test_suite = testEvaluate_Prior_LogNormal
  initTestSuite;

%
% Lognormal distribution:
%
% Parameters:
%  mu
%  sigma > 0.0
%
% Support:
%  0.0 < x < +Inf
%
% http://en.wikipedia.org/wiki/Lognormal_distribution
%

% Test standard functionality
function testStandard
  Order = [ 0, 1, 2 ];
  Expected = [ 0.39894228, -1.0, 0.0 ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_LogNormal([ 0.0, 1.0 ], 1.0, Order(i)));
  end

% Test boundaries of support
function testZeroX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, NaN, NaN ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_LogNormal([ 0.0, 1.0 ], 0.0, Order(i)));
  end

function testInfX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, NaN, NaN ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_LogNormal([ 0.0, 1.0 ], Inf, Order(i)));
  end

% Test parameter checking
%  function testZeroSigma
%    Order = [ 0, 1, 2 ];
%    for i = 1:3
%      assertExceptionThrown(@() Evaluate_Prior_LogNormal([ 0.0, 0.0 ], 1.0, Order), 'LogNormal:Domain_Error');
%    end
%
%  function testOrder
%    assertExceptionThrown(@() Evaluate_Prior_LogNormal([ 1.0, 1.0 ], 1.0, 3), 'Prior:Domain_Error');

% Test Inf propagation
function testInfMu
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, Inf, -Inf ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_LogNormal([ Inf, 1.0 ], 1.0, Order(i)));
  end

function testInfSigma
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, -1.0, 1.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_LogNormal([ 0.0, Inf ], 1.0, Order(i)));
  end

% Test NaN propagation
function testNaNX
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(NaN, Evaluate_Prior_LogNormal([ 0.0, 1.0 ], NaN, Order(i)));
  end

function testNaNMu
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(NaN, Evaluate_Prior_LogNormal([ NaN, 1.0 ], 1.0, Order(i)));
  end

function testNaNSigma
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(NaN, Evaluate_Prior_LogNormal([ 0.0, NaN ], 1.0, Order(i)));
  end
