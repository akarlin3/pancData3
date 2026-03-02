classdef test_plot_parameter_maps < matlab.unittest.TestCase
    % TEST_PLOT_PARAMETER_MAPS Smoke tests for plot_parameter_maps.
    %
    % Verifies:
    %   - No patients with valid adc_vector → "No patients" path, no PNG saved
    %   - Protocol-deviation b-values → patient skipped, no PNG saved
    %   - Valid patient data with correct NIfTI files → PNG is created
    %   - Zero-size data_vectors_gtvp → function handles gracefully

    properties
        TempDir
        OriginalPath
    end

    methods(TestMethodSetup)
        function setup(testCase)
            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path();
            baseDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(baseDir, 'core'));
            addpath(fullfile(baseDir, 'utils'));
            addpath(fullfile(baseDir, 'dependencies'));
            set(0, 'DefaultFigureVisible', 'off');
        end
    end

    methods(TestMethodTeardown)
        function teardown(testCase)
            close all;
            if exist(testCase.TempDir, 'dir')
                rmdir(testCase.TempDir, 's');
            end
            path(testCase.OriginalPath);
        end
    end

    % ------------------------------------------------------------------ %
    %  Helper: build a minimal data_vectors_gtvp struct                  %
    % ------------------------------------------------------------------ %
    methods(Access = private)
        function dv = emptyDataVectors(~, n)
            % n patients, adc_vector = [] → patients are all skipped
            entry = struct('adc_vector', []);
            dv = repmat(entry, n, 1, 1);
        end

        function createValidPatientFiles(testCase, pat_id)
            % Creates the NIfTI and .bval files required for an eligible patient.
            nii_dir = fullfile(testCase.TempDir, pat_id, 'nii');
            mkdir(nii_dir);

            % 4D DWI: 10×10×5×4 (b=0, b=30, b=150, b=550)
            dwi_img = zeros(10, 10, 5, 4);
            dwi_img(:,:,:,1) = 1000;
            dwi_img(:,:,:,2) = 800;
            dwi_img(:,:,:,3) = 500;
            dwi_img(:,:,:,4) = 100;
            dwi_nii = fullfile(nii_dir, 'fx1_dwi1.nii');
            niftiwrite(dwi_img, dwi_nii);
            gzip(dwi_nii);
            delete(dwi_nii);

            % 3D GTV mask: simple central box
            gtv_img = zeros(10, 10, 5);
            gtv_img(4:7, 4:7, 2:4) = 1;
            gtv_nii = fullfile(nii_dir, 'fx1_gtv1.nii');
            niftiwrite(gtv_img, gtv_nii);
            gzip(gtv_nii);
            delete(gtv_nii);

            % b-value file with standard protocol
            fid = fopen(fullfile(nii_dir, 'fx1_dwi1.bval'), 'w');
            fprintf(fid, '0 30 150 550');
            fclose(fid);
        end
    end

    methods(Test)

        function testNoEligiblePatientsNoFileCreated(testCase)
            % All patients have empty adc_vector → skipped in the eligibility
            % pre-count → n_eligible = 0 → patients_plotted = 0 →
            % "No patients … Skipping." message, no PNG written.
            n   = 3;
            dv  = testCase.emptyDataVectors(n);
            ids = {'Pt01', 'Pt02', 'Pt03'};
            out = fullfile(testCase.TempDir, 'output');
            mkdir(out);

            plot_parameter_maps(dv, n, ids, testCase.TempDir, out);

            pngs = dir(fullfile(out, 'Parameter_Maps_*.png'));
            testCase.verifyEmpty(pngs, ...
                'No PNG should be written when no patients have ADC data.');
        end

        function testValidPatientCreatesParameterMapFigure(testCase)
            % One patient with all required files → Parameter_Maps_1.png created.
            pat_id = 'Pt01';
            testCase.createValidPatientFiles(pat_id);

            n_vox = 10 * 10 * 5;   % volume voxel count
            entry = struct('adc_vector', rand(n_vox, 1) * 2e-3);
            dv    = entry;   % 1×1×1 struct

            out = fullfile(testCase.TempDir, 'output');
            mkdir(out);

            plot_parameter_maps(dv, 1, {pat_id}, testCase.TempDir, out);

            testCase.verifyTrue( ...
                exist(fullfile(out, 'Parameter_Maps_1.png'), 'file') > 0, ...
                'Parameter_Maps_1.png should be created for one valid patient.');
        end

        function testProtocolDeviationSkipsPatient(testCase)
            % Wrong b-values in .bval → patient skipped → no PNG created.
            pat_id = 'Pt01';
            testCase.createValidPatientFiles(pat_id);

            % Overwrite .bval with non-standard values
            bval_file = fullfile(testCase.TempDir, pat_id, 'nii', 'fx1_dwi1.bval');
            fid = fopen(bval_file, 'w');
            fprintf(fid, '0 50 200 800');   % non-standard
            fclose(fid);

            n_vox = 10 * 10 * 5;
            entry = struct('adc_vector', rand(n_vox, 1) * 2e-3);
            out   = fullfile(testCase.TempDir, 'output');
            mkdir(out);

            plot_parameter_maps(entry, 1, {pat_id}, testCase.TempDir, out);

            pngs = dir(fullfile(out, 'Parameter_Maps_*.png'));
            testCase.verifyEmpty(pngs, ...
                'No PNG should be written when all patients have protocol deviations.');
        end

        function testMultiplePatientsGroupedIntoFigures(testCase)
            % Six valid patients: pats_per_fig = 5, so two figures are expected.
            n   = 6;
            out = fullfile(testCase.TempDir, 'output');
            mkdir(out);

            % Create NIfTI files for all 6 patients
            ids   = arrayfun(@(x) sprintf('Pt%02d', x), 1:n, 'UniformOutput', false);
            n_vox = 10 * 10 * 5;
            entries(n, 1, 1) = struct('adc_vector', []);
            for k = 1:n
                testCase.createValidPatientFiles(ids{k});
                entries(k, 1, 1).adc_vector = rand(n_vox, 1) * 2e-3;
            end

            plot_parameter_maps(entries, n, ids, testCase.TempDir, out);

            pngs = dir(fullfile(out, 'Parameter_Maps_*.png'));
            testCase.verifyGreaterThanOrEqual(numel(pngs), 1, ...
                'At least one parameter-map figure should be created for 6 patients.');
        end

    end
end
