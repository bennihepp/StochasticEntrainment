function score = simulate_average_and_compute_entrainment_scores(input_period, input_amplitude, options)

    addpath('../');

    S = options;

    [TT, output] = NFkB_TNF_Run(S.Ntrials, S.t0, S.tf, S.dt, S.recordStep, S.volume, ...
        S.input_offset, input_amplitude, input_period);

    %% cutoff transient
    offset_time = S.to;
    offset = find(TT >= offset_time, 1);
    output = output(offset:end, :);
    %% substract mean
    output = output - mean(output, 1);

    %% Fourier spectrum analysis
    for m=1:S.Ntrials
        [omega1, y1] = compute_normalized_fft_truncated(output(:,m)', S.recordStep, 2*pi*S.min_frequency, 2*pi*S.max_frequency);
        if m == 1
            omega = omega1;
            y = zeros(S.Ntrials, length(y1));
        end
        y(m, :) = y1;
    end

    mean_y = mean(y, 1);
    score = compute_entrainment_score(omega, mean_y, input_period, options);

end
