classdef test_scale_td_panel < matlab.unittest.TestCase
    % TEST_SCALE_TD_PANEL Test suite for scale_td_panel function
    %
    % This test verifies the time-dependent panel scaling logic, ensuring:
    % 1. Statistics (mu, sigma) are computed strictly from training rows.
    % 2. Scaling is applied to all rows (train and test) per week.
    % 3. Duplicate rows for the same patient/week are handled correctly (first one used).
    % 4. Special handling for derivative features at week 0.
    % 5. Zero variance and single sample edge cases.

    properties
        OriginalPath
    end

    methods(TestMethodSetup)
        function addUtilsToPath(testCase)
            testCase.OriginalPath = path;
            % Add utils directory to path. Assuming tests/ is one level below root.
            utilsPath = fullfile(pwd, 'utils');
            if ~exist(utilsPath, 'dir')
                % If running from tests/ directory
                utilsPath = fullfile(pwd, '..', 'utils');
            end
            addpath(utilsPath);
        end
    end

    methods(TestMethodTeardown)
        function restorePath(testCase)
            path(testCase.OriginalPath);
        end
    end

    methods(Test)
        function testBasicScaling(testCase)
            % Scenario:
            % Week 1.
            % Train Patients: 1, 2. Values: 10, 20.
            % Test Patient: 3. Value: 30.
            % Expected Mu = (10+20)/2 = 15.
            % Expected Sigma = std([10, 20]) = 7.0711.
            % Scaled Train 1: (10-15)/7.0711 = -0.7071
            % Scaled Train 2: (20-15)/7.0711 = 0.7071
            % Scaled Test 3: (30-15)/7.0711 = 2.1213

            feat_names = {'ADC_Mean'};
            pat_id_td = [1; 2; 3];
            t_start_td = [7; 7; 7]; % All week 1 (ceil(7/7)=1)
            train_pat_ids = [1; 2];

            X_td_raw = [10; 20; 30];

            X_td_scaled = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids);

            mu = 15;
            sig = std([10, 20]);

            expected_vals = (X_td_raw - mu) / sig;

            testCase.verifyEqual(X_td_scaled, expected_vals, 'AbsTol', 1e-4, ...
                'Basic scaling failed for single week scenario');
        end

        function testDerivativeAtWeek0(testCase)
            % Scenario:
            % Week 0 (Day 0).
            % Feature: DeltaADC (contains 'Delta').
            % Logic: Derivative features at Week 0 should have mu=0, sigma=1.
            % This means output = input.

            feat_names = {'DeltaADC'};
            pat_id_td = [1; 2];
            t_start_td = [0; 0];
            train_pat_ids = [1; 2];

            X_td_raw = [100; 200];

            X_td_scaled = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids);

            % Expect no scaling
            testCase.verifyEqual(X_td_scaled, X_td_raw, ...
                'Derivative features at Week 0 should not be scaled (mu=0, sigma=1)');
        end

        function testPatientDuplicateHandling(testCase)
            % Scenario:
            % Week 1.
            % Train Patient 1 has two rows (values 10 and 100).
            % Train Patient 2 has one row (value 20).
            % Logic: Only first row per patient per week is used for stats.
            % Stats calculated from [10, 20]. Value 100 is ignored.
            % Mu = 15, Sigma = 7.0711.
            % All rows are scaled, including the duplicate.

            feat_names = {'ADC_Mean'};
            pat_id_td = [1; 1; 2];
            t_start_td = [7; 7; 7];
            train_pat_ids = [1; 2];

            X_td_raw = [10; 100; 20];

            X_td_scaled = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids);

            mu = 15;
            sig = std([10, 20]);

            expected_vals = (X_td_raw - mu) / sig;

            testCase.verifyEqual(X_td_scaled, expected_vals, 'AbsTol', 1e-4, ...
                'Duplicate rows should not bias statistics calculation');
        end

        function testZeroVariance(testCase)
            % Scenario:
            % Week 1.
            % Train Patients: 1, 2. Values: 10, 10.
            % Mu = 10, Sigma = 0.
            % Logic: Sigma should be clamped to 1.
            % Output = (X - 10) / 1.

            feat_names = {'ADC_Mean'};
            pat_id_td = [1; 2; 3];
            t_start_td = [7; 7; 7];
            train_pat_ids = [1; 2];

            X_td_raw = [10; 10; 50];

            X_td_scaled = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids);

            expected_vals = [0; 0; 40]; % (50-10)/1 = 40

            testCase.verifyEqual(X_td_scaled, expected_vals, 'AbsTol', 1e-4, ...
                'Zero variance should be handled by setting sigma=1');
        end

        function testMultipleWeeks(testCase)
            % Scenario:
            % Week 1: Train [10, 20] -> Mu=15, Sig=7.07
            % Week 2: Train [100, 200] -> Mu=150, Sig=70.71
            % Test rows scaled by respective week stats.

            feat_names = {'ADC_Mean'};
            pat_id_td = [1; 2; 1; 2];
            t_start_td = [7; 7; 14; 14]; % Weeks 1 and 2
            train_pat_ids = [1; 2];

            X_td_raw = [10; 20; 100; 200];

            X_td_scaled = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids);

            mu1 = 15; sig1 = std([10, 20]);
            mu2 = 150; sig2 = std([100, 200]);

            expected_vals = [ ...
                (10-mu1)/sig1; ...
                (20-mu1)/sig1; ...
                (100-mu2)/sig2; ...
                (200-mu2)/sig2 ...
            ];

            testCase.verifyEqual(X_td_scaled, expected_vals, 'AbsTol', 1e-4, ...
                'Multiple weeks should be scaled independently');
        end

        function testSingleTrainingSample(testCase)
            % Scenario:
            % Week 1.
            % Only 1 training patient. Value = 10.
            % Logic: Mu = 10, Sigma = 1.

            feat_names = {'ADC_Mean'};
            pat_id_td = [1; 2];
            t_start_td = [7; 7];
            train_pat_ids = [1]; % Only pat 1 is training

            X_td_raw = [10; 50];

            X_td_scaled = scale_td_panel(X_td_raw, feat_names, pat_id_td, t_start_td, train_pat_ids);

            % Train row: (10-10)/1 = 0
            % Test row: (50-10)/1 = 40
            expected_vals = [0; 40];

            testCase.verifyEqual(X_td_scaled, expected_vals, 'AbsTol', 1e-4, ...
                'Single training sample should result in sigma=1');
        end
    end
end
