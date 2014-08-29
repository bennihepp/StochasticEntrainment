function score = compute_entrainment_score(omega, y, input_period, options)

    S = options;

    if isfield(S, 'entrainment_ratios')
        entrainment_ratios = S.entrainment_ratios;
    elseif isfield(S, 'entrainment_ratio')
        entrainment_ratios = S.entrainment_ratio;
    else
        entrainment_ratios = 1;
    end

    scores = zeros(length(entrainment_ratios), 1);
    for i=1:length(entrainment_ratios)
        ratio = entrainment_ratios(i);

        om_natural = 2 * pi / S.natural_period;
        om_input = 2 * pi / input_period / ratio;
        dom = S.FREQUENCY_NEIGHBOURHOOD_FACTOR * om_natural;

        if abs(om_natural - om_input) < 2 * dom
            score = inf;
        else
            power_total = sum(abs(y).^2);
            power_input = compute_spectrum_power(omega, y, om_input, dom);
            power_input_harmonics = 0;
            for n=2:S.MAX_HARMONIC_N
                power_input_harmonics = power_input_harmonics + compute_spectrum_power(omega, y, om_input * n, dom);
            end
            if power_input >= S.MIN_HARMONICS_POWER_THRESHOLD * power_input_harmonics
                power_input = power_input + power_input_harmonics;
            end

            score = power_input / power_total;

        end

        scores(i) = score;

    end

    score = max(scores);

end
