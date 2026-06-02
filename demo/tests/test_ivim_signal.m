function test_ivim_signal()
% TEST_IVIM_SIGNAL  Analytic checks on the IVIM forward model.
    b = [0; 30; 150; 550];
    S0 = 1000; D = 1.1e-3; f = 0.12; Dstar = 20e-3;

    % At b=0 the model must return exactly S0 (both exponentials are 1).
    S = ivim_signal(b, S0, D, f, Dstar);
    assert(abs(S(1) - S0) < 1e-9, 'S(b=0) must equal S0');

    % Output shape matches the b shape.
    assert(isequal(size(S), size(b)), 'output must match b shape');

    % Strictly decreasing in b for positive D, D* and 0<f<1.
    assert(all(diff(S) < 0), 'signal must decay monotonically with b');

    % f=0 reduces to a pure mono-exponential exp(-b*D).
    S_f0 = ivim_signal(b, S0, D, 0, Dstar);
    assert(max(abs(S_f0 - S0*exp(-b*D))) < 1e-9, 'f=0 must give mono-exp in D');

    % f=1 reduces to exp(-b*D*).
    S_f1 = ivim_signal(b, S0, D, 1, Dstar);
    assert(max(abs(S_f1 - S0*exp(-b*Dstar))) < 1e-9, 'f=1 must give mono-exp in D*');

    % Convex combination: S = f*S0*exp(-bD*) + (1-f)*S0*exp(-bD).
    expect = f*S0*exp(-b*Dstar) + (1-f)*S0*exp(-b*D);
    assert(max(abs(S - expect)) < 1e-9, 'must equal the two-compartment sum');

    % f outside [0,1] is rejected.
    threw = false;
    try; ivim_signal(b, S0, D, 1.5, Dstar); catch; threw = true; end
    assert(threw, 'f>1 must error');
end
