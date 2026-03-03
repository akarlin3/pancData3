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

    if nargin < 2, label = 'Progress'; end
    count = 0;
    callback = @on_complete;

    function on_complete(~)
        count = count + 1;
        text_progress_bar(count, total, label);
    end
end
