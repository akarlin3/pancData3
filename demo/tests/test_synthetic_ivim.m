function test_synthetic_ivim()
% TEST_SYNTHETIC_IVIM  Cohort generator: structure, reproducibility, physics,
% and end-to-end recovery through the REAL pipeline fitter.

    % Small, fast cohort.
    cfg = struct('n_patients', 6, 'n_fractions', 5, 'snr', 60, 'seed', 7);
    c1 = synthetic_ivim(cfg);

    % --- structure & labelling ---
    assert(c1.meta.synthetic == true, 'meta.synthetic must be true');
    nb = numel(c1.meta.bvalues);
    sc = c1.patients(1).fraction(1);
    assert(ndims(sc.dwi) == 4 && size(sc.dwi,4) == nb, 'dwi must be 4-D (x,y,z,b)');
    assert(islogical(sc.mask) && isequal(size(sc.mask), size(c1.gtv_mask)), 'mask must be logical GTV-shaped');
    assert(c1.meta.n_resistant_voxels > 0, 'cohort must contain a resistant core');

    % --- reproducibility under a fixed seed ---
    c2 = synthetic_ivim(cfg);
    assert(isequal(c1.patients(1).fraction(1).dwi, c2.patients(1).fraction(1).dwi), ...
        'same seed must reproduce identical phantom data');
    c3 = synthetic_ivim(setfield(cfg, 'seed', 8)); %#ok<SFLD>
    assert(~isequal(c1.patients(1).fraction(1).dwi, c3.patients(1).fraction(1).dwi), ...
        'different seed must change the data');

    % --- resistant core has lower true D than the ordinary tumour ---
    reg = sc.region; gm = c1.gtv_mask;
    D_core    = mean(sc.D_true(reg == 2 & gm));
    D_ordinary= mean(sc.D_true(reg == 1 & gm));
    assert(D_core < D_ordinary, 'resistant core must have lower true D (restricted diffusion)');

    % --- longitudinal response: a responder's mean true D rises Fx1 -> FxN,
    %     a non-responder's stays ~flat ---
    rp = find(arrayfun(@(x) x.responder, c1.patients), 1);
    nr = find(arrayfun(@(x) ~x.responder, c1.patients), 1);
    if ~isempty(rp)
        d1 = mean(c1.patients(rp).fraction(1).D_true(gm));
        dN = mean(c1.patients(rp).fraction(cfg.n_fractions).D_true(gm));
        assert(dN > d1, 'responder true D must rise across fractions');
    end
    if ~isempty(nr)
        d1 = mean(c1.patients(nr).fraction(1).D_true(gm));
        dN = mean(c1.patients(nr).fraction(cfg.n_fractions).D_true(gm));
        assert(abs(dN - d1) < 1e-4, 'non-responder true D must stay ~flat');
    end

    % --- end-to-end: the REAL fitter recovers D within tolerance at high SNR ---
    rec = fit_synthetic_scan(sc.dwi, c1.meta.bvalues, sc.mask, 100);
    D_truth_mean = mean(sc.D_true(gm));
    D_rec_mean   = nanmean(rec.D);
    rel = abs(D_rec_mean - D_truth_mean) / D_truth_mean;
    assert(rel < 0.20, sprintf('recovered mean D off by %.1f%% at SNR=60 (expect <20%%)', 100*rel));
end
