function c = cvpartition(n, varargin)
    % Mock for cvpartition
    c = struct();
    c.NumTestSets = 1;
    if ischar(n) || isstring(n)
        % e.g. cvpartition(group, 'Leaveout')
        groups = varargin{1};
        c.NumTestSets = length(unique(groups));
        % Basic stub
    else
        % N, 'KFold', k
        c.NumTestSets = 5;
    end
end
