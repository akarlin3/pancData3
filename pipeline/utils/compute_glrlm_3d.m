function [sre, lre, gln, rln, rp] = compute_glrlm_3d(quantized, mask, n_levels)
%COMPUTE_GLRLM_3D  3D Gray-Level Run-Length Matrix features at 13 directions.
%   Computes SRE, LRE, GLN, RLN, RP averaged over 13 unique 3D directions.

    % 13 unique directions in a 3x3x3 neighbourhood
    directions = [
        0  0  1;   % z
        0  1  0;   % y
        1  0  0;   % x
        0  1  1;   % yz diagonal
        0  1 -1;   % yz anti-diagonal
        1  0  1;   % xz diagonal
        1  0 -1;   % xz anti-diagonal
        1  1  0;   % xy diagonal
        1 -1  0;   % xy anti-diagonal
        1  1  1;   % body diagonal
        1  1 -1;   % body anti-diagonal
        1 -1  1;   % body anti-diagonal
        1 -1 -1;   % body anti-diagonal
    ];

    [nr, nc, ns] = size(quantized);
    max_run = max([nr, nc, ns]);

    R = zeros(n_levels, max_run);

    for d = 1:size(directions, 1)
        dr = directions(d, 1);
        dc = directions(d, 2);
        ds = directions(d, 3);

        visited = false(nr, nc, ns);

        for r = 1:nr
            for c = 1:nc
                for s = 1:ns
                    if visited(r,c,s) || quantized(r,c,s) == 0 || ~mask(r,c,s)
                        continue;
                    end

                    gl = quantized(r, c, s);
                    run_len = 1;
                    visited(r, c, s) = true;

                    cr = r + dr; cc = c + dc; cs = s + ds;
                    while cr >= 1 && cr <= nr && cc >= 1 && cc <= nc && ...
                          cs >= 1 && cs <= ns && mask(cr,cc,cs) && ...
                          quantized(cr,cc,cs) == gl
                        visited(cr, cc, cs) = true;
                        run_len = run_len + 1;
                        cr = cr + dr; cc = cc + dc; cs = cs + ds;
                    end

                    if gl >= 1 && gl <= n_levels && run_len >= 1 && run_len <= max_run
                        R(gl, run_len) = R(gl, run_len) + 1;
                    end
                end
            end
        end
    end

    [sre, lre, gln, rln, rp] = glrlm_features(R, sum(mask(:)));
end
