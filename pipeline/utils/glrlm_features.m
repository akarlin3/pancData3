function [sre, lre, gln, rln, rp] = glrlm_features(R, n_voxels)
%GLRLM_FEATURES  Compute SRE, LRE, GLN, RLN, RP from a run-length matrix.

    total_runs = sum(R(:));
    if total_runs == 0
        sre = NaN; lre = NaN; gln = NaN; rln = NaN; rp = NaN;
        return;
    end

    n_gl = size(R, 1);
    max_run = size(R, 2);
    j = 1:max_run;

    % SRE = sum R(i,j)/j^2 / total_runs
    j2_inv = 1 ./ (j.^2);
    sre = sum(R * j2_inv') / total_runs;

    % LRE = sum j^2 * R(i,j) / total_runs
    j2 = j.^2;
    lre = sum(R * j2') / total_runs;

    % GLN = sum_i (sum_j R(i,j))^2 / total_runs
    gl_sums = sum(R, 2);
    gln = sum(gl_sums.^2) / total_runs;

    % RLN = sum_j (sum_i R(i,j))^2 / total_runs
    rl_sums = sum(R, 1);
    rln = sum(rl_sums.^2) / total_runs;

    % RP = total_runs / n_voxels
    rp = total_runs / max(n_voxels, 1);
end
