classdef test_decision_curve_analysis < matlab.unittest.TestCase
% TEST_DECISION_CURVE_ANALYSIS  Tests for decision curve analysis.
%
%   Verifies: (1) perfect predictor has net benefit >= treat-all,
%   (2) random predictor approximates treat-all, (3) all-same-class
%   returns zero net benefit for model.

    properties
        OriginalPath
        TempDir
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
            set(0, 'DefaultFigureVisible', 'off');
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            close all;
            path(testCase.OriginalPath);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods(Test)
        function testPerfectPredictorBeatsTreatAll(testCase)
            % A perfect predictor should have net benefit >= treat-all
            % at all thresholds (where both are defined).
            rng(42);
            n = 100;
            y_true = [ones(30, 1); zeros(70, 1)];
            y_pred = [ones(30, 1); zeros(70, 1)];  % perfect predictions

            thresholds = 0.1:0.1:0.9;
            results = decision_curve_analysis(y_true, y_pred, thresholds);

            % At all thresholds, model NB should be >= treat-all NB
            for i = 1:numel(thresholds)
                testCase.verifyGreaterThanOrEqual( ...
                    results.net_benefit_model(i), ...
                    results.net_benefit_treat_all(i) - 1e-10, ...
                    sprintf('Perfect predictor should beat treat-all at threshold %.1f', thresholds(i)));
            end
        end

        function testRandomPredictorApproximatesTreatAll(testCase)
            % A random noise predictor should approximate the treat-all line.
            rng(42);
            n = 1000;
            y_true = [ones(300, 1); zeros(700, 1)];
            y_pred = rand(n, 1);  % random predictions

            thresholds = 0.2:0.1:0.8;
            results = decision_curve_analysis(y_true, y_pred, thresholds);

            % Net benefit should be close to treat-all (within tolerance)
            for i = 1:numel(thresholds)
                diff = abs(results.net_benefit_model(i) - results.net_benefit_treat_all(i));
                testCase.verifyLessThan(diff, 0.15, ...
                    'Random predictor NB should approximate treat-all.');
            end
        end

        function testAllSameClassZeroNetBenefit(testCase)
            % When all patients have the same class (no events),
            % model net benefit should be 0 at most thresholds.
            n = 50;
            y_true = zeros(n, 1);  % all non-events
            y_pred = rand(n, 1) * 0.5;

            thresholds = 0.1:0.1:0.9;
            results = decision_curve_analysis(y_true, y_pred, thresholds);

            % With no events, TP=0 always, so NB <= 0
            for i = 1:numel(thresholds)
                testCase.verifyLessThanOrEqual(results.net_benefit_model(i), 0 + 1e-10, ...
                    'All-negative outcome should have NB <= 0.');
            end
        end

        function testTreatNoneIsZero(testCase)
            % Treat-none strategy should always have zero net benefit.
            y_true = [1; 0; 1; 0; 1; 0; 0; 0];
            y_pred = rand(8, 1);
            thresholds = 0:0.1:1;
            results = decision_curve_analysis(y_true, y_pred, thresholds);

            testCase.verifyEqual(results.net_benefit_treat_none, ...
                zeros(1, numel(thresholds)), ...
                'Treat-none net benefit should be zero at all thresholds.');
        end

        function testDCAPlotGenerated(testCase)
            % Verify DCA plot is saved when output_folder is provided.
            y_true = [ones(10,1); zeros(20,1)];
            y_pred = [0.8*ones(10,1); 0.2*ones(20,1)];
            results = decision_curve_analysis(y_true, y_pred, 0:0.1:1, ...
                testCase.TempDir, 'Standard', 'Fx2');

            fig_path = fullfile(testCase.TempDir, 'dca_Standard_Fx2.png');
            testCase.verifyTrue(exist(fig_path, 'file') > 0, ...
                'DCA plot should be saved to output folder.');
        end
    end
end
