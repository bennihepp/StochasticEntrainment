function within_arnold_tongue = is_average_within_arnold_tongue(input_period, input_amplitude, options)

    score = simulate_average_and_compute_entrainment_scores(input_period, input_amplitude, options);
    if isnan(score)
        within_arnold_tongue = true;
    else
        within_arnold_tongue = score >= options.ENTRAINMENT_THRESHOLD;
    end;

end
