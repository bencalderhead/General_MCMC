function test_suite = testEvaluate_Prior_Gamma
  initTestSuite;

%
% Gamma distribution:
%
% Parameters:
%  shape > 0.0
%  scale > 0.0
%
% Support:
%  0.0 < x < +Inf
%
% http://en.wikipedia.org/wiki/Gamma_distribution
%

% Test standard functionality
function testStandard
  Order = [ 0, 1, 2 ];
  Expected = [ 0.36787944, -1.0, 0.0 ];
  for i=1:3
    assertElementsAlmostEqual(Expected(i), Evaluate_Prior_Gamma([ 1.0, 1.0 ], 1.0, Order(i)));
  end

% Test boundaries of support
function testZeroX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, -Inf, -Inf ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Gamma([ 1.0, 1.0 ], 0.0, Order(i)));
  end

function testInfX
  Order = [ 0, 1, 2 ];
  Expected = [ 0.0, -1.0, 0.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Gamma([ 1.0, 1.0 ], Inf, Order(i)));
  end

% Test parameter checking
%  function testZeroShape
%    Order = [ 0, 1, 2 ];
%    for i = 1:3
%      assertExceptionThrown(@() Evaluate_Prior_Gamma([ 0.0, 1.0 ], 1.0, Order), 'Gamma:Domain_Error');
%    end
%
%  function testInfShape
%    Order = [ 0, 1, 2 ];
%    Expected = [ 0.0, Inf, -Inf ];
%    for i=1:3
%      assertEqual(Expected(i), Evaluate_Prior_Gamma([ Inf, 1.0 ], 1.0, Order(i)));
%    end
%
%  function testZeroScale
%    Order = [ 0, 1, 2 ];
%    for i = 1:3
%      assertExceptionThrown(@() Evaluate_Prior_Gamma([ 1.0, 0.0 ], 1.0, Order(i)), 'Gamma:Domain_Error');
%    end
%
%  function testInfScale
%    Order = [ 0, 1, 2 ];
%    Expected = [ 0.0, 0.0, 0.0 ];
%    for i=1:3
%      assertEqual(Expected(i), Evaluate_Prior_Gamma([ 1.0, Inf ], 1.0, Order(i)));
%    end
%
%  function testOrder
%    assertExceptionThrown(@() Evaluate_Prior_Gamma([ 1.0, 1.0 ], 1.0, 3), 'Prior:Domain_Error');

% Test NaN propagation
function testNaNX
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(NaN, Evaluate_Prior_Gamma([ 1.0, 1.0 ], NaN, Order(i)));
  end

function testNaNShape
  Order = [ 0, 1, 2 ];
  for i=1:3
    assertEqual(NaN, Evaluate_Prior_Gamma([ NaN, 1.0 ], 1.0, Order(i)));
  end

function testNaNScale
  Order = [ 0, 1, 2 ];
  Expected = [ NaN, NaN, 0.0 ];
  for i=1:3
    assertEqual(Expected(i), Evaluate_Prior_Gamma([ 1.0, NaN ], 1.0, Order(i)));
  end
