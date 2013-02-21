function test_suite = testNormal_LogPDF
  initTestSuite;

%
% Multivariate Normal distribution in log space:
%
% Parameters:
%  mu(n)
%  sigma(n*n) > 0.0
%
% Support:
%  0.0 < x(n) < +Inf
%
% http://en.wikipedia.org/wiki/Normal_distribution
%

% Test standard functionality
function testStandard1
  X = ones(10, 1);
  Mu = ones(10, 1);
  Sigma = eye(10);
  assertElementsAlmostEqual(-9.1894, Normal_LogPDF(X, Mu, Sigma), 'relative', 0.00002);

function testStandard2
  X = ones(10, 1) * 1.5;
  Mu = ones(10, 1);
  Sigma = eye(10) * 5;
  assertElementsAlmostEqual(-17.4866, Normal_LogPDF(X, Mu, Sigma), 'relative', 0.002);

function testStandard3
  X = ones(10, 1) * 0.0001;
  Mu = ones(10, 1);
  Sigma = eye(10);
  assertElementsAlmostEqual(-14.1884, Normal_LogPDF(X, Mu, Sigma), 'relative', 0.0001);

% Test parameter checking
function testNonPosDefSigma
  X = ones(10, 1);
  Mu = ones(10, 1);
  Sigma = eye(10) * -1.0;
  assertExceptionThrown(@() Normal_LogPDF(X, Mu, Sigma), 'MATLAB:posdef');

function testNonSquareSigma
  X = ones(10, 1);
  Mu = ones(10, 1);
  Sigma = eye(10, 8);
  assertExceptionThrown(@() Normal_LogPDF(X, Mu, Sigma), 'MATLAB:square');

function testSizeMismatchX
  X = ones(8, 1);
  Mu = ones(10, 1);
  Sigma = eye(10);
  assertExceptionThrown(@() Normal_LogPDF(X, Mu, Sigma), 'MATLAB:dimagree');

function testSizeMismatchMu
  X = ones(10, 1);
  Mu = ones(12, 1);
  Sigma = eye(10);
  assertExceptionThrown(@() Normal_LogPDF(X, Mu, Sigma), 'MATLAB:dimagree');

function testSizeMismatchSigma
  X = ones(10, 1);
  Mu = ones(10, 1);
  Sigma = eye(8);
  assertExceptionThrown(@() Normal_LogPDF(X, Mu, Sigma), 'MATLAB:dimagree');

% Test NaN propagation
function testNaNX
  X = ones(10, 1);
  Mu = ones(10, 1);
  Sigma = eye(10);
  X(1) = NaN;
  assertEqual(NaN, Normal_LogPDF(X, Mu, Sigma));

function testNaNMu
  X = ones(10, 1);
  Mu = ones(10, 1);
  Sigma = eye(10);
  Mu(1) = NaN;
  assertEqual(NaN, Normal_LogPDF(X, Mu, Sigma));

function testNaNSigma
  X = ones(10, 1);
  Mu = ones(10, 1);
  Sigma = eye(10);
  Sigma(1, 1) = NaN;
  assertExceptionThrown(@() Normal_LogPDF(X, Mu, Sigma), 'MATLAB:posdef');
