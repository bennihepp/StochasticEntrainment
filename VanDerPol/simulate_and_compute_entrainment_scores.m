function score = simulate_and_compute_entrainment_scores(input_period, input_amplitude, options)

    addpath('../');

    S = options;

%     additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
%     multiplicative_forcing_func = @(t, x) 0;

    [TT, output] = VanDerPol_Run_Java(S.Ntrials, S.t0, S.tf, S.dt, S.recordStep, S.volume, S.input_offset, input_amplitude, input_period);
%     [TT, output] = VanDerPol_Run(S.Ntrials, S.t0, S.tf, S.dt, S.volume, ...
%         additive_forcing_func, multiplicative_forcing_func);

    offset_time = S.to;
    offset = find(TT >= offset_time, 1);
    output = output(offset:end, :);

    %% substract mean
    output = output - repmat(mean(output, 1), [size(output, 1), 1]);

    %% Fourier spectrum analysis
    scores = zeros(S.Ntrials, 1);
    for m=1:S.Ntrials
        [omega, y] = compute_normalized_fft_truncated(output(:,m)', S.recordStep, 2*pi*S.min_frequency, 2*pi*S.max_frequency);
        scores(m) = compute_entrainment_score(omega, y, input_period, options);
    end

    score = max(scores);

end
