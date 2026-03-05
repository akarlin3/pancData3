function callback = parfor_progress(total, label)
%PARFOR_PROGRESS Creates a progress callback for use with parallel.pool.DataQueue.
%
%   callback = parfor_progress(total, label)
%
%   Returns a function handle suitable for afterEach(dq, callback).
%   Each invocation increments an internal counter and displays progress.
%
%   Usage:
%       dq = parallel.pool.DataQueue;
%       afterEach(dq, parfor_progress(N, 'Processing patients'));
%       parfor j = 1:N
%           % ... work ...
%           send(dq, j);
%       end
%
%   Inputs:
%       total - Total number of expected completions
%       label - Descriptive label string
%
%   Outputs:
%       callback - Function handle @(~) that tracks and displays progress
%
% --- Analytical Rationale ---
% Processing a pancreatic DWI cohort involves per-patient DICOM conversion,
% model fitting (ADC + IVIM), and optional deep learning denoising. Each
% patient takes 1-10 minutes depending on scan count and model complexity.
% With 30+ patients, the full pipeline can run for hours. Without progress
% feedback, it is impossible to distinguish a stalled pipeline from one that
% is simply slow.
%
% MATLAB's parfor does not natively support progress tracking because loop
% iterations execute out of order on parallel workers. The DataQueue
% mechanism provides a thread-safe channel for workers to notify the client.
% This function creates a closure that maintains a monotonically increasing
% counter on the client side, converting asynchronous worker completions
% into a sequential progress display.
%
% The closure pattern (nested function capturing 'count') is used instead of
% a persistent variable to allow multiple independent progress bars in the
% same MATLAB session (e.g., one for patient processing, another for DIR
% registration).

    if nargin < 2, label = 'Progress'; end

    % --- Closure State ---
    % 'count' is captured by the nested function and persists across calls.
    % It is incremented each time a parfor worker sends a completion signal.
    % The variable lives on the client (not the workers), so no
    % synchronization is needed -- DataQueue serializes delivery.
    count = 0;
    callback = @on_complete;

    function on_complete(~)
        % The argument from send(dq, j) is ignored (~) because we only need
        % to know that a worker finished, not which iteration completed.
        % Iteration indices arrive out of order in parfor, so tracking the
        % count of completions is more meaningful than tracking which
        % specific patient finished.
        count = count + 1;
        text_progress_bar(count, total, label);
    end
end
