function test_suite = testdlngevm
initTestSuite;

function testStandard
  a = [ -Inf, 0 ];
  B = [ 4
        3 ];
  Expected = 3;
  Actual = dlngevm(a, B);
  assertElementsAlmostEqual(Expected, Actual);
