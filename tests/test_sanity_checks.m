classdef test_sanity_checks < matlab.unittest.TestCase

    properties
        OriginalPath
        OriginalPwd
        TempDir
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            testCase.OriginalPath = path;
            testCase.OriginalPwd = pwd;

            % Add core to path
            coreDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core');
            addpath(coreDir);

            % Create and move to temp dir
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            cd(testCase.TempDir);
        end
    end

    methods(TestMethodTeardown)
        function teardownEnvironment(testCase)
            cd(testCase.OriginalPwd);
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    methods
        function [gtvp, gtvn, summary] = createMockData(~, nPat, nTp)
            % Helper to create consistent mock data

            % Initialize summary_metrics
            summary = struct();
            summary.id_list = cell(nPat, 1);
            summary.mrn_list = cell(nPat, 1);
            for i = 1:nPat
                summary.id_list{i} = sprintf('Pat%02d', i);
                summary.mrn_list{i} = sprintf('MRN%02d', i);
            end

            % Initialize metrics arrays [nPat x nTp x 1]
            summary.adc_mean = rand(nPat, nTp, 1);
            summary.d_mean = rand(nPat, nTp, 1);
            summary.f_mean = rand(nPat, nTp, 1);
            summary.dstar_mean = rand(nPat, nTp, 1);
            summary.d95_gtvp = rand(nPat, nTp); % RT Dose metrics

            % Initialize data vectors
            % Structure with fields: adc_vector, d_vector, f_vector, dstar_vector, dose_vector
            % Size: [nPat x nTp x 1]

            % Create an empty struct with fields first to ensure correct array creation
            s = struct('adc_vector', [], 'd_vector', [], 'f_vector', [], ...
                       'dstar_vector', [], 'dose_vector', []);
            gtvp = repmat(s, nPat, nTp, 1);
            gtvn = repmat(s, nPat, nTp, 1);

            for j = 1:nPat
                for k = 1:nTp
                    nVox = 100; % Arbitrary number of voxels

                    % Fill GTVp
                    gtvp(j,k,1).adc_vector = rand(nVox, 1);
                    gtvp(j,k,1).d_vector = rand(nVox, 1);
                    gtvp(j,k,1).f_vector = rand(nVox, 1);
                    gtvp(j,k,1).dstar_vector = rand(nVox, 1);
                    gtvp(j,k,1).dose_vector = rand(nVox, 1);

                    % Fill GTVn (structure is same, data can be different)
                    gtvn(j,k,1).adc_vector = rand(nVox, 1);
                    gtvn(j,k,1).d_vector = rand(nVox, 1);
                    gtvn(j,k,1).f_vector = rand(nVox, 1);
                    gtvn(j,k,1).dstar_vector = rand(nVox, 1);
                    gtvn(j,k,1).dose_vector = rand(nVox, 1);
                end
            end
        end
    end

    methods(Test)
        function testHappyPath(testCase)
            [gtvp, gtvn, summary] = testCase.createMockData(2, 3);

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary);

            testCase.verifyTrue(is_valid, 'Happy path should be valid');
            testCase.verifyTrue(contains(msg, 'Passed'), 'Message should indicate success');
        end

        function testConvergenceWarnings(testCase)
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % Introduce Inf
            gtvp(1,1,1).adc_vector(1) = Inf;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary);

            % Convergence issues are warnings, not failures
            testCase.verifyTrue(is_valid, 'Convergence warnings should not invalidate run');
            testCase.verifyTrue(contains(msg, 'convergence warnings'), 'Message should count warnings');
        end

        function testAlignmentMismatch(testCase)
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % Make dose vector different length
            gtvp(1,1,1).dose_vector(end+1) = 0.5;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary);

            testCase.verifyFalse(is_valid, 'Alignment mismatch should invalidate run');
            testCase.verifyTrue(contains(msg, 'Failed'), 'Message should indicate failure');
            testCase.verifyTrue(contains(msg, 'alignment checks'), 'Message should mention alignment');
        end

        function testDoseNaNFailure(testCase)
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % Make dose vector mostly NaN (>10%)
            nVox = numel(gtvp(1,1,1).dose_vector);
            gtvp(1,1,1).dose_vector(1:floor(0.2*nVox)) = NaN;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary);

            testCase.verifyFalse(is_valid, 'Excessive NaN in dose should invalidate run');
            testCase.verifyTrue(contains(msg, 'Failed'), 'Message should indicate failure');
        end
    end
end
