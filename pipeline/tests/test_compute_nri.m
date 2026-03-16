classdef test_compute_nri < matlab.unittest.TestCase
% TEST_COMPUTE_NRI  Tests for NRI, cNRI, and IDI computation.
%
%   Verifies: (1) perfect reclassification yields NRI near 2.0,
%   (2) identical models yield NRI=0 and IDI=0, (3) edge cases.

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addPaths(testCase)
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'utils'));
        end
    end

    methods(TestMethodTeardown)
        function restorePaths(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)
        function testPerfectReclassificationNRI(testCase)
            % A model that perfectly reclassifies should have NRI near 2.0.
            % Events move up (NRI_events=1), non-events move down (NRI_nonevents=1).
            n = 50;
            y_true = [ones(25, 1); zeros(25, 1)];

            % Old model: everyone in medium risk (0.3)
            prob_old = 0.3 * ones(n, 1);

            % New model: events get high prob, non-events get low prob
            prob_new = [0.9 * ones(25, 1); 0.05 * ones(25, 1)];

            % Categories: [0, 0.2, 0.5, 1.0]
            results = compute_nri(y_true, prob_old, prob_new);

            testCase.verifyGreaterThan(results.nri, 1.5, ...
                'Perfect reclassification should have NRI close to 2.0.');
            testCase.verifyEqual(results.nri_events, 1.0, 'AbsTol', 1e-10, ...
                'All events should reclassify upward.');
            testCase.verifyEqual(results.nri_nonevents, 1.0, 'AbsTol', 1e-10, ...
                'All non-events should reclassify downward.');
        end

        function testIdenticalModelsZeroNRI(testCase)
            % Identical models should produce NRI=0 and IDI=0.
            rng(42);
            n = 40;
            y_true = [ones(15, 1); zeros(25, 1)];
            probs = rand(n, 1);

            results = compute_nri(y_true, probs, probs);

            testCase.verifyEqual(results.nri, 0, 'AbsTol', 1e-10, ...
                'Identical models should have NRI=0.');
            testCase.verifyEqual(results.idi, 0, 'AbsTol', 1e-10, ...
                'Identical models should have IDI=0.');
            testCase.verifyEqual(results.cnri, 0, 'AbsTol', 1e-10, ...
                'Identical models should have cNRI=0.');
        end

        function testEdgeCaseNoEvents(testCase)
            % When there are no events, NRI should still compute without error.
            n = 20;
            y_true = zeros(n, 1);
            prob_old = rand(n, 1) * 0.5;
            prob_new = rand(n, 1) * 0.5;

            results = compute_nri(y_true, prob_old, prob_new);

            testCase.verifyTrue(isfinite(results.nri) || results.nri == 0, ...
                'NRI should be computable with no events.');
        end

        function testOutputFields(testCase)
            % Verify all expected output fields are present.
            y_true = [1; 0; 1; 0; 0];
            prob_old = [0.6; 0.3; 0.7; 0.2; 0.4];
            prob_new = [0.8; 0.1; 0.9; 0.15; 0.3];

            results = compute_nri(y_true, prob_old, prob_new);

            testCase.verifyTrue(isfield(results, 'nri'), 'Missing nri field.');
            testCase.verifyTrue(isfield(results, 'nri_events'), 'Missing nri_events field.');
            testCase.verifyTrue(isfield(results, 'nri_nonevents'), 'Missing nri_nonevents field.');
            testCase.verifyTrue(isfield(results, 'nri_p'), 'Missing nri_p field.');
            testCase.verifyTrue(isfield(results, 'cnri'), 'Missing cnri field.');
            testCase.verifyTrue(isfield(results, 'cnri_ci'), 'Missing cnri_ci field.');
            testCase.verifyTrue(isfield(results, 'idi'), 'Missing idi field.');
            testCase.verifyTrue(isfield(results, 'idi_ci'), 'Missing idi_ci field.');
        end
    end
end
