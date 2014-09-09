function NFkB_TNF_ArnoldTongue_BinarySearch(output_folder, volume, input_periods, input_amplitude_tolerance, Ntrials, population_average)    
    natural_period = 2.1013;
    ENTRAINMENT_THRESHOLD = 0.9;
    MIN_HARMONICS_POWER_THRESHOLD = 0.0;
    FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
    % MAX_HARMONIC_N = 15;
    MAX_HARMONIC_N = double(intmax());
    % entrainment_ratios = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    entrainment_ratios = 1:2;

    dt = 0.0001;
    recordStep = 100 * dt;

    if volume == inf
        Ntrials = 1;
%     else
%         Ntrials = 100;
    end

    t0 = 0;
    tf = 1000;
    to = (tf - t0) / 5;

    input_offset = 1.0;

    % min_input_amplitude = 0.0;
    % max_input_amplitude = 1.0;
    % input_amplitude_tolerance = 1e-2;
    % input_periods = 0.5:0.05:5.0;

    min_input_amplitude = 0.0;
    max_input_amplitude = 0.1;
%     input_amplitude_tolerance = 1e-3;
%     input_periods = 1.9:0.01:2.8;

    min_frequency = 0.0;
    max_frequency = 50.0;


    S = struct();
    S.volume = volume;
    S.population_average = population_average;
    S.Ntrials = Ntrials;
    S.t0 = t0;
    S.tf = tf;
    S.to = to;
    S.dt = dt;
    S.recordStep = recordStep;
    S.input_offset = input_offset;
    S.natural_period = natural_period;
    S.min_frequency = min_frequency;
    S.max_frequency = max_frequency;
    S.input_periods = input_periods;
    S.min_input_amplitude = min_input_amplitude;
    S.max_input_amplitude = max_input_amplitude;
    S.input_amplitude_tolerance = input_amplitude_tolerance;

    S.ENTRAINMENT_THRESHOLD = ENTRAINMENT_THRESHOLD;
    S.MAX_HARMONIC_N = MAX_HARMONIC_N;
    S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
    S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
    S.entrainment_ratios = entrainment_ratios;


    par_arnold_tongue_borders = zeros(length(S.input_periods), 1);
    par_score = zeros(length(S.input_periods), 1);
    par_scores = zeros(length(S.input_periods), 3);
    par_score_mean = zeros(length(S.input_periods), 1);
    par_score_std = zeros(length(S.input_periods), 1);


%     for i=1:length(input_periods)
    parfor i=1:length(input_periods)
        display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
        input_period = input_periods(i);

        display(['input_period=', num2str(input_period)]);

        lower_amplitude = S.min_input_amplitude;
        upper_amplitude = S.max_input_amplitude;

        upper_amp_score = simulate_and_compute_all_entrainment_scores_(input_period, upper_amplitude, S.population_average, S);
        upper_amp_within_at = is_within_arnold_tongue__(upper_amp_score, S);

        display(['upper_amplitude=', num2str(upper_amplitude), ', score=', num2str(upper_amp_score), ', within_at=', num2str(upper_amp_within_at)]);

        lower_amp_score = simulate_and_compute_all_entrainment_scores_(input_period, lower_amplitude, S.population_average, S);
        lower_amp_within_at = is_within_arnold_tongue__(lower_amp_score, S);

        display(['lower_amplitude=', num2str(lower_amplitude), ', score=', num2str(lower_amp_score), ', within_at=', num2str(lower_amp_within_at)]);

        score = inf;
        scores = inf * ones(3, 1);

        if lower_amp_within_at
            display('lower_amp_within_at');
            arnold_tongue_border = lower_amplitude;
            score = lower_amp_score;
        elseif ~upper_amp_within_at
            display('~upper_amp_within_at');
            arnold_tongue_border = inf;
            score = upper_amp_score;
        else

            while (upper_amplitude - lower_amplitude) >= S.input_amplitude_tolerance
                middle_amplitude = (lower_amplitude + upper_amplitude) / 2.0;
                display(['upper_amp=', num2str(upper_amplitude), ', lower_amp=', num2str(lower_amplitude), ', middle_amp=', num2str(middle_amplitude)]);
                score = simulate_and_compute_all_entrainment_scores_(input_period, middle_amplitude, S.population_average, S);
                middle_amp_within_at = is_within_arnold_tongue__(score, S);
                display(['  score=', num2str(score), ', within_at=', num2str(middle_amp_within_at)]);
                if middle_amp_within_at
                    upper_amplitude = middle_amplitude;
                else
                    lower_amplitude = middle_amplitude;
                end
            end

            display(['middle_amp_within_at=', num2str(middle_amp_within_at)]);
            if middle_amp_within_at
                arnold_tongue_border = middle_amplitude;
            else
                arnold_tongue_border = upper_amplitude;
            end

            display(['arnold_tongue_border=', num2str(arnold_tongue_border)]);
            input_amplitude = arnold_tongue_border;
            for j=1:3
                scores(j) = simulate_and_compute_all_entrainment_scores_(input_period, input_amplitude, S.population_average, S);
                display([' j=', int2str(j), ', score=', num2str(scores(j))]);
            end

        end

        score_mean = mean(scores, 1);
        score_std = std(scores, 1);
        display(['score_mean=', num2str(score_mean), ', score_std=', num2str(score_std)]);

        display(['arnold tongue border at ', num2str(arnold_tongue_border)]);

        par_arnold_tongue_borders(i) = arnold_tongue_border;
        par_score(i) = score;
        par_scores(i, :) = scores;
        par_score_mean(i) = score_mean;
        par_score_std(i) = score_std;
    end

    S.natural_period = natural_period;
    S.input_periods = input_periods;
    S.min_input_amplitude = min_input_amplitude;
    S.max_input_amplitude = max_input_amplitude;
    S.input_amplitude_tolerance = input_amplitude_tolerance;
    S.arnold_tongue_borders = par_arnold_tongue_borders;
    S.score = par_score;
    S.scores = par_scores;
    S.score_mean = par_score_mean;
    S.score_std = par_score_std;

    date_string = datestr(clock());
    filename = [output_folder, 'NFkB_TNF_ArnoldTongue_BinarySearch_volume=', num2str(S.volume), '_population=', num2str(S.population_average), '_', date_string, '.mat'];
    save(filename, '-struct', 'S');
end

function within_arnold_tongue = is_within_arnold_tongue__(score, options)
    if any(isnan(score))
        within_arnold_tongue = true;
    else
%         if length(score) > 1
%             ratio_within_arnold_tongue = sum(scores >= options.ENTRAINMENT_THRESHOLD) / length(scores);
%             within_arnold_tongue = ratio_within_arnold_tongue >= options.MINIMUM_ENTRAINMENT_RATIO;
%         end
        within_arnold_tongue = score >= options.ENTRAINMENT_THRESHOLD;
    end;
end

function score = simulate_and_compute_all_entrainment_scores_(input_period, input_amplitude, population_average, options)
    if population_average
        score = simulate_average_and_compute_entrainment_scores(input_period, input_amplitude, options);
    else
        score = simulate_and_compute_all_entrainment_scores(input_period, input_amplitude, options);
        if length(score) > 1
            score = sort(score);
            score = mean(score(round(length(score) / 2):end));
        end
    end
end
