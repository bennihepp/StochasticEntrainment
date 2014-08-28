function score = simulate_and_compute_entrainment_scores(input_period, input_amplitude, options)

    addpath('../');

    S = options;

    [TT, output] = CircadianClock_Run(S.Ntrials, S.t0, S.tf, S.dt, S.recordStep, S.volume, ...
        S.input_offset, input_amplitude, input_period);

    offset_time = S.to;
    offset = find(TT >= offset_time, 1);
    output = output(offset:end, :);

    %% Fourier spectrum analysis
    scores = zeros(S.Ntrials, 1);
    for m=1:S.Ntrials
        [omega, y] = compute_normalized_fft_truncated(output(:,m)', S.recordStep, 2*pi*S.min_frequency, 2*pi*S.max_frequency);
        scores(m) = compute_entrainment_score(omega, y, input_period, options);
    end

    score = max(scores);

end
