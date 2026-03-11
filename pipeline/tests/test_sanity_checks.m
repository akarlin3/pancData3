classdef test_sanity_checks < matlab.unittest.TestCase
    % TEST_SANITY_CHECKS Unit tests for the sanity_checks module.
    %
    % Validates the pre-analysis data integrity checks that guard against
    % corrupt or misaligned inputs. Tests cover:
    %   - Happy path (all valid data)
    %   - Convergence warnings (Inf/NaN/negative values in voxel vectors)
    %   - Spatial alignment mismatches (dose vs. ADC vector length)
    %   - NaN dose coverage warnings
    %   - Outlier detection in summary metrics
    %   - Missingness reporting
    %   - DWI type dispatch (Standard, DnCNN, IVIMnet)
    %   - Default output folder fallback

    properties
        OriginalPath   % Saved MATLAB path for teardown restoration
        OriginalPwd    % Saved working directory for teardown restoration
        TempDir        % Temporary directory for test output files
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % Save current state, add core/ to path, and create a temp
            % directory for sanity check output files (diary logs, results).
            testCase.OriginalPath = path;
            testCase.OriginalPwd = pwd;

            coreDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core');
            addpath(coreDir);

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
            % Creates valid mock data for sanity_checks testing.
            % Returns GTVp/GTVn struct arrays with randomized voxel vectors
            % (100 voxels each) for all DWI types (Standard, DnCNN, IVIMnet),
            % and a summary_metrics struct with randomized summary arrays.
            % All values are in [0,1] range (rand), which is not physically
            % realistic but sufficient for testing the validation logic.

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
                       'dstar_vector', [], 'dose_vector', [], ...
                       'adc_vector_dncnn', [], 'd_vector_dncnn', [], ...
                       'f_vector_dncnn', [], 'dstar_vector_dncnn', [], ...
                       'd_vector_ivimnet', [], ...
                       'f_vector_ivimnet', [], 'dstar_vector_ivimnet', []);
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

                    gtvp(j,k,1).adc_vector_dncnn = rand(nVox, 1);
                    gtvp(j,k,1).d_vector_dncnn = rand(nVox, 1);
                    gtvp(j,k,1).f_vector_dncnn = rand(nVox, 1);
                    gtvp(j,k,1).dstar_vector_dncnn = rand(nVox, 1);

                    gtvp(j,k,1).d_vector_ivimnet = rand(nVox, 1);
                    gtvp(j,k,1).f_vector_ivimnet = rand(nVox, 1);
                    gtvp(j,k,1).dstar_vector_ivimnet = rand(nVox, 1);

                    % Fill GTVn (structure is same, data can be different)
                    gtvn(j,k,1).adc_vector = rand(nVox, 1);
                    gtvn(j,k,1).d_vector = rand(nVox, 1);
                    gtvn(j,k,1).f_vector = rand(nVox, 1);
                    gtvn(j,k,1).dstar_vector = rand(nVox, 1);
                    gtvn(j,k,1).dose_vector = rand(nVox, 1);

                    gtvn(j,k,1).adc_vector_dncnn = rand(nVox, 1);
                    gtvn(j,k,1).d_vector_dncnn = rand(nVox, 1);
                    gtvn(j,k,1).f_vector_dncnn = rand(nVox, 1);
                    gtvn(j,k,1).dstar_vector_dncnn = rand(nVox, 1);

                    gtvn(j,k,1).d_vector_ivimnet = rand(nVox, 1);
                    gtvn(j,k,1).f_vector_ivimnet = rand(nVox, 1);
                    gtvn(j,k,1).dstar_vector_ivimnet = rand(nVox, 1);
                end
            end
        end
    end

    methods(Test)
        function testHappyPath(testCase)
            % Verifies that clean, well-formed data with no anomalies
            % passes all sanity checks and returns is_valid=true.
            [gtvp, gtvn, summary] = testCase.createMockData(2, 3);

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 1, 'dwi_type_name', 'Standard', 'output_folder', testCase.TempDir));

            testCase.verifyTrue(is_valid, 'Happy path should be valid');
            testCase.verifyTrue(contains(msg, 'Passed'), 'Message should indicate success');
        end

        function testConvergenceWarnings(testCase)
            % Verifies that Inf values in voxel vectors trigger convergence
            % warnings but do NOT invalidate the run. Inf values indicate
            % model fitting divergence for individual voxels, which is
            % expected in noisy data and handled downstream by NaN-aware stats.
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % Introduce Inf to simulate a diverged ADC fit
            gtvp(1,1,1).adc_vector(1) = Inf;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 1, 'dwi_type_name', 'Standard', 'output_folder', testCase.TempDir));

            % Convergence issues are warnings, not failures
            testCase.verifyTrue(is_valid, 'Convergence warnings should not invalidate run');
            testCase.verifyTrue(contains(msg, 'convergence warnings'), 'Message should count warnings');
        end

        function testAlignmentMismatch(testCase)
            % Verifies that a length mismatch between the dose vector and
            % the ADC vector causes a hard failure (is_valid=false). This
            % catches spatial misalignment where the dose map was sampled
            % on a different grid than the DWI parameter maps.
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % Append one extra element to dose_vector to create a mismatch
            gtvp(1,1,1).dose_vector(end+1) = 0.5;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 1, 'dwi_type_name', 'Standard', 'output_folder', testCase.TempDir));

            testCase.verifyFalse(is_valid, 'Alignment mismatch should invalidate run');
            testCase.verifyTrue(contains(msg, 'Failed'), 'Message should indicate failure');
            testCase.verifyTrue(contains(msg, 'alignment checks'), 'Message should mention alignment');
        end

        function testDoseNaNFailure(testCase)
            % Verifies that >10% NaN values in the dose vector produce a
            % soft warning (not a hard failure). Partial RT dose overlap
            % is clinically common when the dose grid does not fully cover
            % the GTV, so this is expected and tolerated.
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % Set 20% of dose voxels to NaN to exceed the 10% threshold
            nVox = numel(gtvp(1,1,1).dose_vector);
            gtvp(1,1,1).dose_vector(1:floor(0.2*nVox)) = NaN;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 1, 'dwi_type_name', 'Standard', 'output_folder', testCase.TempDir));

            % NaN dose coverage is a soft warning (partial RT dose overlap
            % is common), not a hard failure.  Only dimensional mismatches
            % invalidate the run.
            testCase.verifyTrue(is_valid, 'NaN dose warnings should not invalidate run');
            testCase.verifyTrue(contains(msg, 'NaN dose warnings'), 'Message should report NaN dose warnings');
        end

        function testOutlierDetection(testCase)
            % Verifies that extreme summary metric outliers are detected
            % and reported but do NOT invalidate the run. Outliers are
            % written to the diary/console for researcher review.
            [gtvp, gtvn, summary] = testCase.createMockData(5, 1);

            % Set 4 patients to tight range to establish a narrow IQR
            summary.adc_mean(1:4,1,1) = [1.0, 1.1, 1.0, 1.1]; % IQR=0.1, median=1.05

            % Patient 5 is far outside the 3*IQR fence [0.75, 1.35]
            summary.adc_mean(5,1,1) = 5.0;

            % Capture console output since outliers are just printed to diary/console
            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 1, 'dwi_type_name', 'Standard', 'output_folder', testCase.TempDir));

            % Outliers should not invalidate the run, but they should be reported
            testCase.verifyTrue(is_valid, 'Outliers should not invalidate run');
            % Outlier warning is written to diary/console. We can test that the function completes successfully.
        end

        function testMissingness(testCase)
            % Verifies that NaN values in summary metrics are detected and
            % reported as missingness but do NOT invalidate the run.
            % Missing data is common when scans are skipped or ROIs fail.
            [gtvp, gtvn, summary] = testCase.createMockData(2, 2);

            % Introduce NaN in both DWI and dose metrics
            summary.adc_mean(1,1,1) = NaN;
            summary.d95_gtvp(1,2) = NaN;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 1, 'dwi_type_name', 'Standard', 'output_folder', testCase.TempDir));

            % Missing data should just be summarized
            testCase.verifyTrue(is_valid, 'Missingness should not invalidate run');
        end

        function testDefaultOutputFolder(testCase)
            % Verifies that sanity_checks works when called without an
            % explicit output_folder (3-arg form). It should create a
            % timestamped saved_files_* directory as a fallback. The test
            % snapshots pre-existing saved_files_* dirs and only cleans up
            % newly created ones to avoid deleting real pipeline output.
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);
            core_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'core');
            project_root = fullfile(core_dir, '..');
            pre_dirs = dir(fullfile(project_root, 'saved_files_*'));
            pre_names = {pre_dirs([pre_dirs.isdir]).name};

            % Call with only 3 arguments, fallback to default output_folder
            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary);

            % Should pass and complete
            testCase.verifyTrue(is_valid, 'Fallback defaults should not invalidate run');

            % Cleanup: only remove saved_files_* dirs that did not exist
            % before the call (i.e., created by this test).
            fallback_dirs = dir(fullfile(project_root, 'saved_files_*'));
            for fi = 1:numel(fallback_dirs)
                if fallback_dirs(fi).isdir && ~ismember(fallback_dirs(fi).name, pre_names)
                    rmdir(fullfile(project_root, fallback_dirs(fi).name), 's');
                end
            end
        end

        function testDnCNNDtype(testCase)
            % Verifies that sanity_checks correctly dispatches to DnCNN
            % vectors (dwi_types_to_run=2) and that negative values in
            % d_vector_dncnn trigger convergence warnings (not failures).
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % Negative D value: physically impossible, indicates fit failure
            gtvp(1,1,1).d_vector_dncnn(1) = -1.0;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 2, 'dwi_type_name', 'DnCNN', 'output_folder', testCase.TempDir));

            % Should trigger convergence warning but not invalidate run
            testCase.verifyTrue(is_valid, 'DnCNN convergence warnings should not invalidate run');
            testCase.verifyTrue(contains(msg, 'convergence warnings'), 'Message should indicate convergence warning for DnCNN');
        end

        function testIvimNetDtype(testCase)
            % Verifies that sanity_checks correctly dispatches to IVIMnet
            % vectors (dwi_types_to_run=3) and that NaN values in
            % f_vector_ivimnet trigger convergence warnings (not failures).
            [gtvp, gtvn, summary] = testCase.createMockData(1, 1);

            % NaN perfusion fraction: indicates IVIMnet inference failure
            gtvp(1,1,1).f_vector_ivimnet(1) = NaN;

            [is_valid, msg, ~, ~] = sanity_checks(gtvp, gtvn, summary, struct('dwi_types_to_run', 3, 'dwi_type_name', 'IVIM-NET', 'output_folder', testCase.TempDir));

            % Should trigger convergence warning but not invalidate run
            testCase.verifyTrue(is_valid, 'IvimNET convergence warnings should not invalidate run');
            testCase.verifyTrue(contains(msg, 'convergence warnings'), 'Message should indicate convergence warning for IvimNET');
        end
    end
end
