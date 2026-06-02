function [sre, lre, gln, rln, rp] = compute_glrlm(quantized, mask, n_levels)
%COMPUTE_GLRLM  2D Gray-Level Run-Length Matrix features at 4 directions.
%   Computes SRE, LRE, GLN, RLN, RP averaged over 4 in-plane directions
%   (0, 45, 90, 135 degrees).

    directions = [0 1; 1 0; 1 1; 1 -1];
    [nrows, ncols] = size(quantized);
    max_run = max(nrows, ncols);

    % Accumulate GLRLM across all directions
    R = zeros(n_levels, max_run);

    for d = 1:size(directions, 1)
        dr = directions(d, 1);
        dc = directions(d, 2);

        visited = false(nrows, ncols);

        for r = 1:nrows
            for c = 1:ncols
                if visited(r, c) || quantized(r, c) == 0 || ~mask(r, c)
                    continue;
                end

                gl = quantized(r, c);
                run_len = 1;
                visited(r, c) = true;

                % Extend run in the current direction
                cr = r + dr;
                cc = c + dc;
                while cr >= 1 && cr <= nrows && cc >= 1 && cc <= ncols && ...
                      mask(cr, cc) && quantized(cr, cc) == gl
                    visited(cr, cc) = true;
                    run_len = run_len + 1;
                    cr = cr + dr;
                    cc = cc + dc;
                end

                if gl >= 1 && gl <= n_levels && run_len >= 1 && run_len <= max_run
                    R(gl, run_len) = R(gl, run_len) + 1;
                end
            end
        end
    end

    [sre, lre, gln, rln, rp] = glrlm_features(R, sum(mask(:)));
end
