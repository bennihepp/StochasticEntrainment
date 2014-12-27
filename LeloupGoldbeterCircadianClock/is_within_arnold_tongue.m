function within_arnold_tongue = is_within_arnold_tongue(input_period, input_amplitude, options)

    scores = simulate_and_compute_all_entrainment_scores(input_period, input_amplitude, options);
    score = max(scores);
    if isnan(score)
        within_arnold_tongue = true;
    else
        within_arnold_tongue = score >= options.ENTRAINMENT_THRESHOLD;
    end;

end
