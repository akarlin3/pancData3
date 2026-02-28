classdef test_dvh < matlab.unittest.TestCase
    % TEST_DVH
    %
    % Verify the functionality of the dose-volume histogram (DVH) calculation.

    methods(TestMethodSetup)
        function setupPaths(testCase)
            % Add necessary paths
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(baseDir, 'dependencies')));
        end
    end

    methods(Test)

        function testBasicCalculation(testCase)
            % Create a simple 3D uniform dose map
            % A 10x10x10 cube
            dose_map = zeros(10, 10, 10);

            % Set dose in a specific region
            % Region 1: 5x5x5 cube with dose 50 (125 voxels)
            % Region 2: 5x5x5 cube with dose 100 (125 voxels)
            dose_map(1:5, 1:5, 1:5) = 50;
            dose_map(6:10, 6:10, 6:10) = 100;

            % Create a structure mask that covers both regions
            struct_mask = false(10, 10, 10);
            struct_mask(1:5, 1:5, 1:5) = true;
            struct_mask(6:10, 6:10, 6:10) = true;
            % Total mask volume is 250 voxels

            % Voxel dimensions in cm
            dims = [0.1, 0.1, 0.1]; % 1 voxel = 0.001 cc
            % Total mask volume = 250 * 0.001 = 0.25 cc

            % Calculate DVH
            [dvhparams, dvh_values] = dvh(dose_map, struct_mask, dims, 201, ...
                'Dperc', 90, ...
                'Dvol', 0.1, ...
                'Vperc', 60, ...
                'Vvol', 40, ...
                'Normalize', true);

            % Verify table columns dynamically generated via dot syntax
            % Vperc for 60 Gy: only the 100 Gy region receives >= 60. That's 125 voxels, which is 50% of the volume.
            testCase.verifyEqual(dvhparams.('V60Gy (%)'), 50, 'AbsTol', 1);

            % Vvol for 40 Gy: both regions receive >= 40. That's 250 voxels = 0.25 cc.
            testCase.verifyEqual(dvhparams.('V40.0Gy (cc)'), 0.25, 'AbsTol', 0.01);

            % For dvh_values, size should be [201, 2]
            testCase.verifyEqual(size(dvh_values), [201, 2]);

            % Max dose in dvh_values(:,1) should be max(dose_map) = 100
            testCase.verifyEqual(max(dvh_values(:, 1)), 100);

            % The first bin dose is 0, volume should be 100%
            testCase.verifyEqual(dvh_values(1, 1), 0);
            testCase.verifyEqual(dvh_values(1, 2), 100, 'AbsTol', 1e-4);
        end

        function testNormalizationFalse(testCase)
            dose_map = zeros(10, 10, 10);
            dose_map(1:5, 1:5, 1:5) = 50;
            struct_mask = false(10, 10, 10);
            struct_mask(1:5, 1:5, 1:5) = true;
            dims = [0.1, 0.1, 0.1];

            [~, dvh_values] = dvh(dose_map, struct_mask, dims, 100, 'Normalize', false);

            % Max volume should be 125 * 0.001 = 0.125 cc
            % first bin is dose 0, so 100% of struct receives > 0, which is 0.125
            testCase.verifyEqual(dvh_values(1, 2), 0.125, 'AbsTol', 0.001);
        end

        function testEmptyMask(testCase)
            dose_map = ones(5, 5, 5) * 50;
            struct_mask = false(5, 5, 5);
            dims = [0.1, 0.1, 0.1];

            [~, dvh_values] = dvh(dose_map, struct_mask, dims, 100, 'Normalize', true);

            % Should handle without crashing, dvh_values(:, 2) should be NaNs
            testCase.verifyTrue(all(isnan(dvh_values(:, 2))));
        end
    end
end