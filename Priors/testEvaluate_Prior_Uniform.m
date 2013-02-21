function test_suite = testEvaluate_Prior_Uniform
  initTestSuite;

%
% Uniform distribution:
%
% Parameters:
%  -Inf < a < b < Inf
%
% Support:
%  a <= x <= b
%
% http://en.wikipedia.org/wiki/Uniform_distribution_(continuous)
%

% Test standard functionality
function testStandard
  Order = [ 0, 1, 2 ];
  Expected = [ 1.0, 0.0, 0.0 ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_Uniform([ 0.0, 1.0 ], 0.5, Order(i)));
  end

% Test boundaries of support
function testUpperX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, -Inf, -Inf ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Uniform([ 0.0, 1.0 ], 1.001, Order(i)));
  end

function testLowerX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, -Inf, -Inf ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Uniform([ 0.0, 1.0 ], -0.001, Order(i)));
  end

% Test parameter checking
%  function testUpperLower
%    Order = [ 0, 1, 2 ];
%    for i = 1:3
%      assertExceptionThrown(@() Evaluate_Prior_Uniform([ 1.0, 0.0 ], 0.0, Order), 'Uniform:Domain_Error');
%    end
%
%  function testUpperLowerEq
%    Order = [ 0, 1, 2 ];
%    for i = 1:3
%      assertExceptionThrown(@() Evaluate_Prior_Uniform([ 0.0, 0.0 ], 0.0, Order), 'Uniform:Domain_Error');
%    end
%
%  function testOrder
%    assertExceptionThrown(@() Evaluate_Prior_Uniform([ 0.0, 1.0 ], 1.0, 4), 'Prior:Domain_Error');

% Test Inf propagation
function testInf
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(0.0, Evaluate_Prior_Uniform([ 0.0, Inf ], 0.0, Order(i)));
  end

% Test NaN propagation
function testNaNX
  Order = [ 0, 1, 2 ];
  Expected = [ 1.0, 0.0, 0.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Uniform([ 0.0, 1.0 ], NaN, Order(i)));
  end

function testNaNMu
  Order = [ 0, 1, 2 ];
  Expected = [ NaN, 0.0, 0.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Uniform([ NaN, 1.0 ], 0.0, Order(i)));
  end

function testNaNSigma
  Order = [ 0, 1, 2 ];
  Expected = [ NaN, 0.0, 0.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Uniform([ 0.0, NaN ], 0.0, Order(i)));
  end
