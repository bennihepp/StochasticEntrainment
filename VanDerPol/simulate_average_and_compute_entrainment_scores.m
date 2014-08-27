function score = simulate_average_and_compute_entrainment_scores(input_period, input_amplitude, options)

    addpath('../');

    S = options;

    additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
    multiplicative_forcing_func = @(t, x) 0;

    [TT, output] = VanDerPol_Run(S.Ntrials, S.t0, S.tf, S.dt, S.volume, ...
        additive_forcing_func, multiplicative_forcing_func);

    offset_time = S.to;
    offset = find(TT >= offset_time, 1);
    output = output(offset:end, :);

    %% Fourier spectrum analysis
    for m=1:S.Ntrials
        [omega1, y1] = compute_normalized_fft_truncated(output(:,m)', S.dt, 2*pi*S.min_frequency, 2*pi*S.max_frequency);
        if m == 1
            omega = omega1;
            y = zeros(S.Ntrials, length(y1));
        end
        y(m, :) = y1;
    end

    mean_y = mean(y, 1);
    score = compute_entrainment_score(omega, mean_y, input_period, options);

end
