addpath('core');
addpath('tests');
testCase = test_compute_summary_metrics();
testCase.createMockInputs();
testCase.testBasicCalculation();
testCase.testSubvolumeCalculation();
disp('Done');
