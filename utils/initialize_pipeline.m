function [resolved_config_path, tests_passed, tests_timestamp] = initialize_pipeline( ...
    pipeline_dir, config_path, steps_to_run, master_output_folder, ...
    tests_passed_in, tests_timestamp_in)
% INITIALIZE_PIPELINE  Pre-flight initialization for the DWI pipeline.
%
%   [resolved_config_path, tests_passed, tests_timestamp] = initialize_pipeline( ...
%       pipeline_dir, config_path, steps_to_run, master_output_folder, ...
%       tests_passed_in, tests_timestamp_in)
%
%   Performs all initialization tasks that precede the main pipeline loop:
%     1. Adds core/, utils/, dependencies/ to the MATLAB path
%     2. Resolves relative config paths against pipeline_dir
%     3. Runs the pre-flight test suite (once per session, with staleness check)
%     4. Verifies required MATLAB toolbox licenses
%
%   Inputs:
%     pipeline_dir         - Absolute path to the pancData3 repository root
%     config_path          - Path to config.json (absolute or relative to pipeline_dir)
%     steps_to_run         - Cell array of pipeline step names
%     master_output_folder - Parent output folder (may be empty)
%     tests_passed_in      - Current value of the persistent tests_passed flag
%                            (empty or logical)
%     tests_timestamp_in   - Current value of the persistent tests_passed_timestamp
%                            (empty or datenum)
%
%   Outputs:
%     resolved_config_path - Fully resolved path to config.json
%     tests_passed         - Updated tests_passed flag
%     tests_timestamp      - Updated tests_passed_timestamp
%
%   Errors:
%     Throws 'PipelineAborted:TestFailure' if the test suite fails.
%     Throws 'InitializationError:MissingToolbox' if required toolboxes are absent.

    % --- 1) Add folders to the MATLAB path ---
    addpath(fullfile(pipeline_dir, 'core'));
    addpath(fullfile(pipeline_dir, 'utils'));
    % Octave compatibility shims provide reimplementations of MATLAB-specific
    % functions (e.g., niftiread, cvpartition, nanmean) so the pipeline can
    % run on GNU Octave for sites without MATLAB licenses. These shims are
    % only loaded when running under Octave to avoid shadowing native MATLAB.
    if exist('OCTAVE_VERSION', 'builtin')
        addpath(fullfile(pipeline_dir, 'utils', 'octave_compat'));
    end
    % The dependencies/ folder contains third-party IVIM fitting algorithms
    % (segmented, Bayesian), ADC fitting, DnCNN denoising, and dose-volume
    % histogram tools. These are treated as read-only external libraries.
    addpath(fullfile(pipeline_dir, 'dependencies'));

    % --- 2) Resolve config path ---
    % If the config path is relative (e.g., 'config.json'), try resolving
    % it against the pipeline root directory. This supports both absolute
    % paths and calls from within the repository directory.
    resolved_config_path = config_path;
    if ~isfile(resolved_config_path) && isfile(fullfile(pipeline_dir, resolved_config_path))
        resolved_config_path = fullfile(pipeline_dir, resolved_config_path);
    end

    % --- 3) Pre-flight test suite ---
    tests_passed = tests_passed_in;
    tests_timestamp = tests_timestamp_in;

    % Check for pre-flight skip via environment variable (used by CI or
    % automated scripts that have already validated the codebase) or via
    % the config.json skip_tests flag (used during development).
    skip_preflight = strcmp(getenv('SKIP_PIPELINE_PREFLIGHT'), '1');
    if ~skip_preflight
        try
            pf_raw = fileread(resolved_config_path);
            pf_cfg = jsondecode(pf_raw);
            if isfield(pf_cfg, 'skip_tests') && pf_cfg.skip_tests
                skip_preflight = true;
            end
        catch
            % If config can't be read here, let parse_config handle the error later
        end
    end
    if ismember('test', steps_to_run)
        % Invalidate cached test result if any test file was modified since
        % the last successful run (supports interactive development).
        if tests_passed && ~isempty(tests_timestamp)
            test_files_info = dir(fullfile(pipeline_dir, 'tests', '**', '*.m'));
            if ~isempty(test_files_info)
                latest_mod = max([test_files_info.datenum]);
                if latest_mod > tests_timestamp
                    tests_passed = false;
                    fprintf('  💡 Test files modified since last pre-flight; re-running.\n');
                end
            end
        end
        if ~skip_preflight && (isempty(tests_passed) || ~tests_passed)
            try
                fprintf('⚙️ [Pre-flight] Running unit tests before pipeline...\n');
                % Run tests directly instead of via run_all_tests.m to:
                %  - Exclude integration tests (test_dwi_pipeline, test_modularity)
                %    that call run_dwi_pipeline themselves, which would be circular.
                %  - Skip the CodeCoveragePlugin (pre-flight only needs pass/fail).
                %  - Keep pre-flight fast (~seconds vs. minutes for full coverage).
                % The excluded integration tests are run separately by
                % execute_all_workflows via run_all_tests.m before any pipeline
                % calls, so they are still exercised in full workflow runs.
                tests_dir = fullfile(pipeline_dir, 'tests');
                addpath(tests_dir);
                pf_suite = matlab.unittest.TestSuite.fromFolder(tests_dir, 'IncludingSubfolders', true);
                % Exclude integration tests that invoke the pipeline
                integration_names = {'test_dwi_pipeline', 'test_modularity'};
                keep = true(size(pf_suite));
                for ii = 1:numel(pf_suite)
                    for jj = 1:numel(integration_names)
                        if startsWith(pf_suite(ii).Name, integration_names{jj})
                            keep(ii) = false;
                        end
                    end
                end
                pf_suite = pf_suite(keep);
                % Capture test output to its own diary file
                if ~isempty(master_output_folder) && exist(master_output_folder, 'dir')
                    pf_diary = fullfile(master_output_folder, 'preflight_tests_output.log');
                    diary(pf_diary);
                end
                % Suppress figure rendering during tests to avoid hundreds of
                % plot windows spawning on headless cluster nodes or during
                % batch execution. Tests validate plot generation logic
                % without needing visible output.
                old_fig_vis = get(0, 'DefaultFigureVisible');
                set(0, 'DefaultFigureVisible', 'off');
                pf_results = run(pf_suite);
                set(0, 'DefaultFigureVisible', old_fig_vis);
                if ~isempty(master_output_folder) && exist(master_output_folder, 'dir')
                    diary off;
                end
                if any([pf_results.Failed])
                    error('PreFlight:TestFailure', '%d test(s) failed.', sum([pf_results.Failed]));
                end
                tests_passed = true;
                tests_timestamp = now;
                fprintf('      ✅ %d unit tests passed.\n', numel(pf_results));
            catch ME
                if ~isempty(master_output_folder) && exist(master_output_folder, 'dir')
                    diary off;
                end
                tests_passed = false;
                fprintf('❌ Test suite failed: %s\n', ME.message);
                error('PipelineAborted:TestFailure', ...
                    'Pipeline aborted because the test suite did not pass.');
            end
        else
            fprintf('⏭️ [Pre-flight] Test suite already passed this session or configured to skip. Skipping.\n');
        end
    else
        fprintf('⏭️ [Pre-flight] Skipping test step.\n');
    end

    % --- 4) Toolbox license checks ---
    if ~exist('OCTAVE_VERSION', 'builtin')
        if ~license('test', 'Statistics_Toolbox')
            error('InitializationError:MissingToolbox', ...
                'The "Statistics and Machine Learning Toolbox" is required but not installed or licensed.');
        end

        if ~license('test', 'Image_Toolbox')
            error('InitializationError:MissingToolbox', ...
                'The "Image Processing Toolbox" is required but not installed or licensed.');
        end
    end
end
