function filepath = discover_gtv_file(folder, patterns, index)
    if ischar(patterns) || isstring(patterns)
        patterns = {char(patterns)};
    end
    filepath = '';
    gtv_search = dir(fullfile(folder, '__nonexistent__'));  % empty struct array with dir fields
    single_gtv_search = gtv_search;
    for p = 1:length(patterns)
        pat = patterns{p};
        tmp = dir(fullfile(folder, [pat int2str(index) '*.mat']));
        if ~isempty(tmp); gtv_search = [gtv_search; tmp]; end
        tmp = dir(fullfile(folder, [pat '.mat']));
        if ~isempty(tmp); single_gtv_search = [single_gtv_search; tmp]; end
    end
    if isempty(gtv_search)
        % fallback: sometimes only 1 mask exists for all repeats
        if ~isempty(single_gtv_search)
            filepath = fullfile(single_gtv_search(1).folder, single_gtv_search(1).name);
        end
    else
        filepath = fullfile(gtv_search(1).folder, gtv_search(1).name);
    end
end