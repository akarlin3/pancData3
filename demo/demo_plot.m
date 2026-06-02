function demo_plot(spec, outpath)
% DEMO_PLOT  Render a labelled multi-series 2-D plot to a PNG.
%
% A small, dependency-light plotting helper for the synthetic demo. It tries
% NATIVE MATLAB/Octave graphics first (the normal path in a real MATLAB
% install); if that fails — e.g. a headless CI box or an Octave build with no
% working GL/text renderer — it falls back to writing a gnuplot script and
% shelling out to `gnuplot`. Either way the demo produces the same figure, so
% the headline figures can be regenerated anywhere.
%
% spec fields:
%   .title, .xlabel, .ylabel   strings
%   .series                    struct array, each with:
%        .x, .y                vectors (same length)
%        .yerr   (optional)    symmetric error-bar half-heights
%        .label                legend entry
%        .kind   (optional)    'line' | 'points' | 'linepoints'  (default 'line')
%   .identity (optional)       true -> draw a y=x reference line (recovery plots)
%   .caption  (optional)       a one-line "SYNTHETIC" stamp drawn on the figure
%
% outpath: destination .png path.

    if ~isfield(spec, 'caption'); spec.caption = 'SYNTHETIC PHANTOM DATA — not clinical'; end
    if ~isfield(spec, 'identity'); spec.identity = false; end

    % Decide the renderer ONCE per session: a real MATLAB/Octave install
    % renders natively; a headless box with a broken GL/text renderer falls
    % back to gnuplot. Probing once avoids repeatedly attempting (and noisily
    % failing) native rendering for every figure.
    if native_renderer_ok()
        if try_native(spec, outpath); return; end
    end
    try_gnuplot(spec, outpath);
end

% -------------------------------------------------------------------------
function tf = native_renderer_ok()
    persistent cached
    if ~isempty(cached); tf = cached; return; end
    tf = false; f = [];
    try
        f = figure('visible', 'off');
        axes('parent', f);            % throws on a broken headless renderer
        tf = true;
    catch
        tf = false;
    end
    if ~isempty(f) && ishandle(f); try; close(f); catch; end; end
    cached = tf;
end

% -------------------------------------------------------------------------
function ok = try_native(spec, outpath)
% Native MATLAB/Octave rendering. Wrapped so a broken headless renderer
% (which throws at axes creation) cleanly hands off to the gnuplot fallback.
    ok = false;
    f = [];
    try
        f = figure('visible', 'off', 'position', [100 100 760 520]);
        ax = axes('parent', f); hold(ax, 'on');
        if spec.identity
            allx = []; for s = spec.series, allx = [allx; s.x(:)]; end %#ok<AGROW>
            ally = []; for s = spec.series, ally = [ally; s.y(:)]; end %#ok<AGROW>
            lo = min([allx; ally]); hi = max([allx; ally]);
            plot(ax, [lo hi], [lo hi], 'k--', 'displayname', 'identity (truth = recovered)');
        end
        for s = spec.series
            kind = 'line'; if isfield(s, 'kind') && ~isempty(s.kind); kind = s.kind; end
            if isfield(s, 'yerr') && ~isempty(s.yerr)
                errorbar(ax, s.x, s.y, s.yerr, 'displayname', s.label, 'linewidth', 1.4);
            else
                switch kind
                    case 'points';     plot(ax, s.x, s.y, 'o', 'displayname', s.label);
                    case 'linepoints'; plot(ax, s.x, s.y, '-o', 'displayname', s.label, 'linewidth', 1.4);
                    otherwise;         plot(ax, s.x, s.y, '-',  'displayname', s.label, 'linewidth', 1.4);
                end
            end
        end
        title(ax, spec.title); xlabel(ax, spec.xlabel); ylabel(ax, spec.ylabel);
        legend(ax, 'location', 'best'); grid(ax, 'on');
        % SYNTHETIC stamp so the figure can never be mistaken for clinical data.
        xl = xlim(ax); yl = ylim(ax);
        text(ax, xl(1) + 0.02*diff(xl), yl(2) - 0.04*diff(yl), spec.caption, ...
            'color', [0.6 0 0], 'fontweight', 'bold', 'fontsize', 9);
        print(f, outpath, '-dpng', '-r110');
        close(f);
        ok = true;
    catch
        if ~isempty(f) && ishandle(f); try; close(f); catch; end; end
        ok = false;
    end
end

% -------------------------------------------------------------------------
function try_gnuplot(spec, outpath)
% gnuplot fallback: write per-series data files + a script, then run gnuplot.
    if isempty(which_exe('gnuplot'))
        warning('demo_plot:noRenderer', ...
            'Native plotting failed and gnuplot is not on PATH. Skipping figure %s.', outpath);
        return;
    end
    tmpdir = [outpath '.gpdir'];
    if ~exist(tmpdir, 'dir'); mkdir(tmpdir); end
    cleaner = onCleanup(@() rmdir_safe(tmpdir));

    plot_cmds = {};
    if spec.identity
        allx = []; ally = [];
        for s = spec.series; allx = [allx; s.x(:)]; ally = [ally; s.y(:)]; end %#ok<AGROW>
        lo = min([allx; ally]); hi = max([allx; ally]);
        idf = fullfile(tmpdir, 'identity.dat');
        write_dat(idf, [lo hi]', [lo hi]', []);
        plot_cmds{end+1} = sprintf('"%s" using 1:2 with lines lc rgb "black" dt 2 title "identity"', idf);
    end
    for i = 1:numel(spec.series)
        s = spec.series(i);
        df = fullfile(tmpdir, sprintf('series%d.dat', i));
        if isfield(s, 'yerr') && ~isempty(s.yerr)
            write_dat(df, s.x(:), s.y(:), s.yerr(:));
            plot_cmds{end+1} = sprintf('"%s" using 1:2:3 with yerrorlines lw 2 title "%s"', df, gp_escape(s.label)); %#ok<AGROW>
        else
            write_dat(df, s.x(:), s.y(:), []);
            kind = 'line'; if isfield(s, 'kind') && ~isempty(s.kind); kind = s.kind; end
            switch kind
                case 'points';     style = 'with points pt 7 ps 0.6';
                case 'linepoints'; style = 'with linespoints lw 2 pt 7';
                otherwise;         style = 'with lines lw 2';
            end
            plot_cmds{end+1} = sprintf('"%s" using 1:2 %s title "%s"', df, style, gp_escape(s.label)); %#ok<AGROW>
        end
    end

    gp = fullfile(tmpdir, 'plot.gp');
    fid = fopen(gp, 'w');
    fprintf(fid, 'set terminal pngcairo size 760,520 font "DejaVu Sans,11"\n');
    fprintf(fid, 'set output "%s"\n', outpath);
    fprintf(fid, 'set title "%s"\n', gp_escape(spec.title));
    fprintf(fid, 'set xlabel "%s"\n', gp_escape(spec.xlabel));
    fprintf(fid, 'set ylabel "%s"\n', gp_escape(spec.ylabel));
    fprintf(fid, 'set grid\nset key outside right top\n');
    fprintf(fid, 'set label "%s" at graph 0.02,0.96 tc rgb "#990000" font ",9"\n', gp_escape(spec.caption));
    fprintf(fid, 'plot %s\n', strjoin(plot_cmds, ', \\\n     '));
    fclose(fid);

    cmd = sprintf('gnuplot "%s"', gp);
    [st, msg] = system(cmd);
    if st ~= 0
        warning('demo_plot:gnuplotFailed', 'gnuplot failed for %s: %s', outpath, msg);
    end
end

% -------------------------------------------------------------------------
function write_dat(path, x, y, e)
    fid = fopen(path, 'w');
    if isempty(e)
        fprintf(fid, '%.10g %.10g\n', [x(:) y(:)]');
    else
        fprintf(fid, '%.10g %.10g %.10g\n', [x(:) y(:) e(:)]');
    end
    fclose(fid);
end

function s = gp_escape(s)
    s = strrep(s, '"', '');
    s = strrep(s, '\', '');
end

function p = which_exe(name)
    [st, out] = system(sprintf('command -v %s 2>/dev/null', name));
    if st == 0; p = strtrim(out); else; p = ''; end
end

function rmdir_safe(d)
    try; if exist(d, 'dir'); rmdir(d, 's'); end; catch; end
end
