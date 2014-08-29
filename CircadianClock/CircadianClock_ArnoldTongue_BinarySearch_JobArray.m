% brutus
% POPULATION_AVERAGE=false
% POPULATION_AVERAGE=true
% bsub -n 1 -R "rusage[mem=1536]" -W 12:00 -J "job[1-81]" -o logs/CircadianClock_ArnoldTongue_BinarySearch_JobArray_%I.out bash CircadianClock_ArnoldTongue_BinarySearch_JobArray.sh "\$LSB_JOBINDEX" output/ $POPULATION_AVERAGE
% bsub -n 1 -R "rusage[mem=1536]" -W 8:00 -J "job_comb" -o logs/CircadianClock_ArnoldTongue_BinarySearch_JobArray_Combine.out bash CircadianClock_ArnoldTongue_BinarySearch_JobArray.sh 0 output/ $POPULATION_AVERAGE
%
% INDEX=71; bsub -n 1 -R "rusage[mem=2048]" -W 16:00 -o logs/CircadianClock_ArnoldTongue_BinarySearch_JobArray_$INDEX.out bash CircadianClock_ArnoldTongue_BinarySearch_JobArray.sh $INDEX output/ $POPULATION_AVERAGE
%

function CircadianClock_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, population_average)

    if nargin < 3
        population_average = false;
    end

    ENTRAINMENT_THRESHOLD = 0.9;
    MAX_HARMONIC_N = 4;
    MIN_HARMONICS_POWER_THRESHOLD = 0.0;
    FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
%     STD_ESTIMATION_SIZE = 3;
    natural_period = 23.7473;

    volume = 1e-20;

    if volume == inf
        Ntrials = 1;
%         Ntrials_levels = [1];
%         Ntrials_std = [0];
        dt = 0.002;
        recordStep = 100 * dt;
    else
        Ntrials = 100;
%         Ntrials_levels = [50, 100, 200, 500, 1000];
%         Ntrials_std = zeros(size(Ntrials_levels));
        dt = 0.002;
        recordStep = 100 * dt;
    end

    disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);

    t0 = 0;
    tf = 200*72;
    to = (tf - t0) / 5;

    input_offset = 1.0;

    min_input_amplitude = 0.0;
    max_input_amplitude = 1.0;
    input_amplitude_tolerance = 1e-2;
    input_periods = 12:0.25:36;

    input_amplitude_tolerance = 2e-1;
    input_periods = 12:1.0:36;

    min_frequency = 0.005;
    max_frequency = 0.5;


    S = struct();
    S.volume = volume;
    S.population_average = population_average;
    S.Ntrials = Ntrials;
%     S.Ntrials_levels = Ntrials_levels;
%     S.Ntrials_std = Ntrials_std;
    S.t0 = t0;
    S.tf = tf;
    S.to = to;
    S.dt = dt;
    S.recordStep = recordStep;
    S.input_offset = input_offset;
    S.natural_period = natural_period;
    S.min_frequency = min_frequency;
    S.max_frequency = max_frequency;
    S.ENTRAINMENT_THRESHOLD = ENTRAINMENT_THRESHOLD;
    S.MAX_HARMONIC_N = MAX_HARMONIC_N;
    S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
    S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
%     S.STD_ESTIMATION_SIZE = STD_ESTIMATION_SIZE;


%     if n == -1
%         %% compute variances for Ntrials_levels
% 
%         input_period = median(input_periods);
%         input_amplitude = 0.5;
%         for level=1:length(Ntrials_levels)
%             S.Ntrials = Ntrials_levels(level);
%             scores = zeros(STD_ESTIMATION_SIZE, 1);
%             for i=1:length(scores)
%                 scores(i) = simulate_and_compute_entrainment_score_(input_period, input_amplitude, population_average, S);
%             end
%             S.Ntrials_std(level) = std(scores);
%         end
% 
%         filename = get_Ntrials_std_filename(filename_prefix);
%         save(filename, '-struct', 'S');
% 
%         return;
%     end
% 
%     filename = get_Ntrials_std_filename(filename_prefix);
%     S = load(filename);

    if n == 0
        %% combine results of all jobs

        arnold_tongue_borders = zeros(length(input_periods), 1);
%         Ntrials = zeros(length(input_periods), 1);
        score = zeros(length(input_periods), 1);
        score_mean = zeros(length(input_periods), 1);
        score_std = zeros(length(input_periods), 1);

        for i=1:length(input_periods)
            display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
            filename = get_filename(filename_prefix, i);
            error_occured = false;
            try
                tmp_S = load(filename);
            catch
                error_occured = true;
                display(['error occured: skipping i=', int2str(i)]);
            end
            if ~error_occured
                arnold_tongue_borders(i) = tmp_S.arnold_tongue_border;
%                 Ntrials(i) = tmp_S.Ntrials;
                score(i) = tmp_S.score;
                score_mean(i) = tmp_S.score_mean;
                score_std(i) = tmp_S.score_std;
            end
        end

        S.natural_period = natural_period;
        S.input_periods = input_periods;
        S.min_input_amplitude = min_input_amplitude;
        S.max_input_amplitude = max_input_amplitude;
        S.input_amplitude_tolerance = input_amplitude_tolerance;
        S.arnold_tongue_borders = arnold_tongue_borders;
        S.score = score;
        S.score_mean = score_mean;
        S.score_std = score_std;
%         S.Ntrials = Ntrials;

        date_string = datestr(clock());
        filename = ['CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=', num2str(volume), '_', date_string, '.mat'];
        save(filename, '-struct', 'S');

%         if ~error_occured
%             for i=1:length(input_periods)
%                 filename = [filename_prefix, 'job_', int2str(n), '.mat'];
%                 display(['deleting file: ', filename]);
% %                 delete(filename);
%             end
%         end

    else
        %% perform a single job

        i = n;

        display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
        input_period = input_periods(i);

        lower_amplitude = min_input_amplitude;
        upper_amplitude = max_input_amplitude;

        upper_amp_within_at = is_within_arnold_tongue_(input_period, upper_amplitude, population_average, S);
        lower_amp_within_at = is_within_arnold_tongue_(input_period, lower_amplitude, population_average, S);

        score = inf;
        score_mean = inf;
        score_std = inf;

        if lower_amp_within_at
            arnold_tongue_border = lower_amplitude;
        elseif ~upper_amp_within_at
            arnold_tongue_border = inf;
        else

            while (upper_amplitude - lower_amplitude) >= input_amplitude_tolerance
                middle_amplitude = (lower_amplitude + upper_amplitude) / 2.0;
                display(['i=', int2str(i), ' of ', int2str(length(input_periods)), ', trying input=', num2str(middle_amplitude)]);
                score = simulate_and_compute_entrainment_score_(input_period, middle_amplitude, population_average, S);
                middle_amp_within_at = is_within_arnold_tongue__(score, S);
                if middle_amp_within_at
                    upper_amplitude = middle_amplitude;
                else
                    lower_amplitude = middle_amplitude;
                end
            end

            if middle_amp_within_at
                arnold_tongue_border = middle_amplitude;
                input_amplitude = middle_amplitude;
            else
                arnold_tongue_border = upper_amplitude;
                input_amplitude = upper_amplitude;
            end

            scores = zeros(3, 1);
            for j=1:3
                scores(j) = simulate_and_compute_entrainment_score_(input_period, input_amplitude, population_average, S);
            end
            score_mean = mean(scores);
            score_std = std(scores);

        end

        display(['i=', int2str(i), ' arnold tongue border at ', num2str(arnold_tongue_border)]);

        S = struct();
        S.arnold_tongue_border = arnold_tongue_border;
        S.score = score;
        S.score_mean = score_mean;
        S.score_std = score_std;
%         S.Ntrials = Ntrials;

        filename = get_filename(filename_prefix, n);
        save(filename, '-struct', 'S');

    end

end

% function filename = get_Ntrials_std_filename(prefix)
%     filename = [prefix, 'Ntrials_std.mat'];
% end

function filename = get_filename(prefix, n)
    filename = [prefix, 'job_', int2str(n), '.mat'];
end

function within_arnold_tongue = is_within_arnold_tongue__(score, options)
    if isnan(score)
        within_arnold_tongue = true;
    else
        within_arnold_tongue = score >= options.ENTRAINMENT_THRESHOLD;
    end;
end

function within_arnold_tongue = is_within_arnold_tongue_(input_period, input_amplitude, population_average, options)
    if population_average
        within_arnold_tongue = is_average_within_arnold_tongue(input_period, input_amplitude, options);
    else
        within_arnold_tongue = is_within_arnold_tongue(input_period, input_amplitude, options);
    end
end

function score = simulate_and_compute_entrainment_score_(input_period, input_amplitude, population_average, options)
    if population_average
        score = simulate_average_and_compute_entrainment_scores(input_period, input_amplitude, options);
    else
        score = simulate_and_compute_entrainment_scores(input_period, input_amplitude, options);
    end
end
