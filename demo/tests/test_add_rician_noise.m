function test_add_rician_noise()
% TEST_ADD_RICIAN_NOISE  Statistical properties of the Rician noise model.
    rng(1234, 'twister');

    % sigma = 0 is the noiseless degenerate case: out = |S|.
    S = [0 100 1000];
    assert(isequal(add_rician_noise(S, 0), abs(S)), 'sigma=0 must return |S|');

    % Output is always non-negative (it is a magnitude).
    sigma = 50;
    out = add_rician_noise(repmat(500, 1, 10000), sigma);
    assert(all(out >= 0), 'Rician output must be non-negative');

    % Noise floor: at true signal S=0 the magnitude is Rayleigh with mean
    % sigma*sqrt(pi/2). Check the empirical mean matches within 3%.
    floor_samples = add_rician_noise(zeros(1, 200000), sigma);
    expected_floor = sigma * sqrt(pi/2);
    rel = abs(mean(floor_samples) - expected_floor) / expected_floor;
    assert(rel < 0.03, sprintf('noise floor mean off by %.1f%%', 100*rel));

    % High SNR: the magnitude mean is approximately the true signal.
    big = 100000;
    hi = add_rician_noise(repmat(big, 1, 50000), sigma);
    assert(abs(mean(hi) - big)/big < 0.001, 'high-SNR mean must approach S');

    % Reproducibility under a fixed seed.
    rng(99, 'twister'); a = add_rician_noise(repmat(300, 1, 100), 20);
    rng(99, 'twister'); c = add_rician_noise(repmat(300, 1, 100), 20);
    assert(isequal(a, c), 'same seed must give identical noise');

    % Negative sigma rejected.
    threw = false;
    try; add_rician_noise(S, -1); catch; threw = true; end
    assert(threw, 'negative sigma must error');
end
