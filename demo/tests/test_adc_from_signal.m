function test_adc_from_signal()
% TEST_ADC_FROM_SIGNAL  ADC estimator round-trip and guards.
    b = [0; 30; 150; 550];
    S0 = 1000;

    % Mono-exponential signal must be recovered nearly exactly (noiseless).
    for adc_true = [0.5e-3, 1.0e-3, 1.8e-3]
        S = S0 * exp(-b * adc_true);
        adc = adc_from_signal(b, S);
        rel = abs(adc - adc_true) / adc_true;
        assert(rel < 1e-3, sprintf('ADC round-trip off by %.3g%% at ADC=%.2g', 100*rel, adc_true));
    end

    % b(1) must be 0.
    threw = false;
    try; adc_from_signal([10;150;550], [1;2;3]); catch; threw = true; end
    assert(threw, 'non-zero first b must error');

    % Non-positive signal returns NaN (log undefined).
    assert(isnan(adc_from_signal(b, [1000;0;100;50])), 'zero signal must give NaN');

    % An increasing signal (non-physical) returns NaN, not a negative ADC.
    assert(isnan(adc_from_signal(b, [100;200;400;800])), 'increasing signal must give NaN');
end
