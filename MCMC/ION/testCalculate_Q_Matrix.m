function test_suite = testCalculate_Q_Matrix
initTestSuite;

function testCastillo_Katz
  ModelName = 'Castillo_Katz';
  RateParas = [ 1, 2, 3, 4 ];
  Expected = [ -1,  1,  0,
                2, -5,  3,
                0,  4, -4 ];
  Actual = Calculate_Q_Matrix(ModelName, RateParas);
  assertVectorsAlmostEqual(Expected, Actual);

function testFive_State
  ModelName = 'Five_State';
  RateParas = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ];
  Expected = [  -8.0,   0.0,     0,   8.0,     0,
                12.0, -22.0,  10.0,     0,     0,
                   0,   9.0, -17.0,   8.0,     0,
                 7.0,     0,   0.0,  -9.0,   2.0,
                   0,     0,     0,   0.0,  -0.0 ];
  Actual = Calculate_Q_Matrix(ModelName, RateParas);
  assertVectorsAlmostEqual(Expected, Actual, 'relative', 1e-7);
