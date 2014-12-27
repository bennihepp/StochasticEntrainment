function within_arnold_tongue = is_within_arnold_tongue(input_period, input_amplitude, options)

    score = simulate_and_compute_entrainment_scores_ODE(input_period, input_amplitude, options);
    if isnan(score)
        within_arnold_tongue = true;
    else
        within_arnold_tongue = score >= options.ENTRAINMENT_THRESHOLD;
    end;

end
