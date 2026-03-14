classdef test_dvh < matlab.unittest.TestCase
    % TEST_DVH Unit tests for the dose-volume histogram (DVH) calculation.
    %
    % Validates dvh.m (in dependencies/) which computes cumulative DVH
    % statistics for a given dose map and structure mask. These metrics are
    % used in the dosimetry pipeline (metrics_dosimetry.m) to characterize
    % radiation dose coverage of the GTV and surrounding anatomy.
    %
    % Tests cover:
    %   - Basic DVH calculation with known dose geometry (V60Gy, V40Gy)
    %   - Absolute volume mode (Normalize=false)
    %   - Edge case: empty structure mask (all-NaN output)
    %
    % Note: dvh.m uses MATLAB's 'arguments' block, which is not supported
    % in Octave, so all tests skip on Octave.

    methods(TestMethodSetup)
        function setupPaths(testCase)
            % Add dependencies/ to path for dvh.m access
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            if exist('OCTAVE_VERSION', 'builtin')
                addpath(fullfile(baseDir, 'dependencies'));
            else
                testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(baseDir, 'dependencies')));
            end
        end
    end

    methods(Test)

        function testBasicCalculation(testCase)
            % Verify DVH statistics for a known two-region dose geometry.
            % Region 1 (125 voxels): uniform 50 Gy
            % Region 2 (125 voxels): uniform 100 Gy
            % Total mask: 250 voxels. This setup makes hand-calculation easy:
            %   V60Gy(%) = 50% (only the 100 Gy region gets >= 60 Gy)
            %   V40Gy(cc) = 0.25 cc (both regions get >= 40 Gy)
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            dose_map = zeros(10, 10, 10);

            % Assign 50 Gy to one octant and 100 Gy to the diagonally opposite octant
            dose_map(1:5, 1:5, 1:5) = 50;     % Region 1: 125 voxels at 50 Gy
            dose_map(6:10, 6:10, 6:10) = 100;  % Region 2: 125 voxels at 100 Gy

            % Structure mask covers both dose regions (250 voxels total)
            struct_mask = false(10, 10, 10);
            struct_mask(1:5, 1:5, 1:5) = true;
            struct_mask(6:10, 6:10, 6:10) = true;
            % Total mask volume is 250 voxels

            % Voxel dimensions in cm (0.1 x 0.1 x 0.1 = 0.001 cc per voxel)
            % Total mask volume = 250 voxels * 0.001 cc = 0.25 cc
            dims = [0.1, 0.1, 0.1];

            % Calculate DVH with 201 bins, requesting multiple DVH metrics:
            %   Dperc=90: D90 (dose covering 90% of volume)
            %   Dvol=0.1: Dose at 0.1 cc absolute volume
            %   Vperc=60: V60Gy (% volume receiving >= 60 Gy)
            %   Vvol=40:  V40Gy (absolute cc receiving >= 40 Gy)
            [dvhparams, dvh_values] = dvh(dose_map, struct_mask, dims, 201, ...
                'Dperc', 90, ...
                'Dvol', 0.1, ...
                'Vperc', 60, ...
                'Vvol', 40, ...
                'Normalize', true);

            % Verify Vperc for 60 Gy: only the 100 Gy region (125/250 = 50%) receives >= 60 Gy
            testCase.verifyEqual(dvhparams.('V60Gy (%)'), 50, 'AbsTol', 1);

            % Verify Vvol for 40 Gy: both regions (250 voxels = 0.25 cc) receive >= 40 Gy
            testCase.verifyEqual(dvhparams.('V40.0Gy (cc)'), 0.25, 'AbsTol', 0.01);

            % dvh_values is [nBins x 2]: column 1 = dose, column 2 = cumulative volume
            testCase.verifyEqual(size(dvh_values), [201, 2]);

            % The dose axis should span from 0 to max(dose_map)=100 Gy
            testCase.verifyEqual(max(dvh_values(:, 1)), 100);

            % At dose=0, 100% of the structure receives at least 0 Gy (cumulative DVH starts at 100%)
            testCase.verifyEqual(dvh_values(1, 1), 0);
            testCase.verifyEqual(dvh_values(1, 2), 100, 'AbsTol', 1e-4);
        end

        function testNormalizationFalse(testCase)
            % Verify that with Normalize=false, dvh_values reports absolute
            % volume in cc rather than percentage. For 125 voxels at 0.001 cc
            % each, the total volume at dose=0 should be 0.125 cc.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
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
            % Edge case: an empty structure mask (no voxels selected) should
            % not crash. The cumulative volume column should be all NaN since
            % there are zero voxels to compute statistics over.
            if exist('OCTAVE_VERSION', 'builtin'); return; end
            dose_map = ones(5, 5, 5) * 50;
            struct_mask = false(5, 5, 5);
            dims = [0.1, 0.1, 0.1];

            [~, dvh_values] = dvh(dose_map, struct_mask, dims, 100, 'Normalize', true);

            % Should handle without crashing, dvh_values(:, 2) should be NaNs
            testCase.verifyTrue(all(isnan(dvh_values(:, 2))));
        end
    end
end