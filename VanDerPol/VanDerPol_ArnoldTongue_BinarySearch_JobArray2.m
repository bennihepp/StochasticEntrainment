% brutus
% POPULATION_AVERAGE=false
% POPULATION_AVERAGE=true
% INPUT_PERIODS=12:0.25:36
% AMPLITUDE_TOLERANCE=1e-2
% INPUT_PERIODS=[24,30,36]
% NTRIALS=100
% AMPLITUDE_TOLERANCE=1e-1
% bsub -n 1 -R "rusage[mem=1536]" -W 8:00 -J "job1a" -o logs/VanDerPol_ArnoldTongue_BinarySearch_JobArray_Combine.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh -1 output/ "$INPUT_PERIODS" $AMPLITUDE_TOLERANCE $NTRIALS $POPULATION_AVERAGE
% bsub -n 1 -R "rusage[mem=2048]" -W 16:00 -w "done(job1a)" -J "job1b[1-3]" -o logs/VanDerPol_ArnoldTongue_BinarySearch_JobArray_%I.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh "\$LSB_JOBINDEX" output/
% bsub -n 1 -R "rusage[mem=1536]" -W 8:00 -w "done(job1b)" -J "job1c" -o logs/VanDerPol_ArnoldTongue_BinarySearch_JobArray_Combine.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh 0 output/
%
% INDEX=71; bsub -n 1 -R "rusage[mem=2048]" -W 16:00 -o logs/VanDerPol_ArnoldTongue_BinarySearch_JobArray_$INDEX.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh $INDEX output/ $POPULATION_AVERAGE
%

function VanDerPol_ArnoldTongue_BinarySearch_JobArray2(n, output_folder, ...
    input_periods, input_amplitude_tolerance, Ntrials, population_average)

%     if n == -2
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
%         filename = get_Ntrials_std_filename(output_folder);
%         save(filename, '-struct', 'S');
% 
%         return;
%     end
% 
%     filename = get_Ntrials_std_filename(output_folder);
%     S = load(filename);

    if n == -1
        ENTRAINMENT_THRESHOLD = 0.9;
        MAX_HARMONIC_N = 4;
        MIN_HARMONICS_POWER_THRESHOLD = 0.0;
        FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
        MINIMUM_ENTRAINMENT_RATIO = 0.5;
    %     STD_ESTIMATION_SIZE = 3;
        natural_period = 23.7473;
        entrainment_ratios = [1, 2];

        volume = 1e-20;

        if volume == inf
    %         Ntrials_levels = [1];
    %         Ntrials_std = [0];
            dt = 0.002;
            recordStep = 100 * dt;
        else
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
        S.input_periods = input_periods;
        S.min_input_amplitude = min_input_amplitude;
        S.max_input_amplitude = max_input_amplitude;
        S.input_amplitude_tolerance = input_amplitude_tolerance;

        S.ENTRAINMENT_THRESHOLD = ENTRAINMENT_THRESHOLD;
        S.MAX_HARMONIC_N = MAX_HARMONIC_N;
        S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
        S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
        S.MINIMUM_ENTRAINMENT_RATIO = MINIMUM_ENTRAINMENT_RATIO;
    %     S.STD_ESTIMATION_SIZE = STD_ESTIMATION_SIZE;
        S.entrainment_ratios = entrainment_ratios;


        %% create MAT-file for all jobs
        info_filename = get_info_filename(output_folder);
        save(info_filename, '-struct', 'S');

        return;
    end

    info_filename = get_info_filename(output_folder);
    S = load(info_filename);

    if n == 0
        %% combine results of all jobs

        arnold_tongue_borders = zeros(length(S.input_periods), 1);
%         Ntrials = zeros(length(input_periods), 1);
        score = zeros(length(S.input_periods), 1);
        score_mean = zeros(length(S.input_periods), 1);
        score_std = zeros(length(S.input_periods), 1);

        for i=1:length(S.input_periods)
            display(['i=', int2str(i), ' out of ', int2str(length(S.input_periods))]);
            filename = get_job_filename(output_folder, i);
            error_occured = false;
            try
                display(['  loading: ', filename]);
                tmp_S = load(filename);
            catch e
                error_occured = true;
                display(['error occured: skipping i=', int2str(i)]);
                display(e);
            end
            if ~error_occured
                arnold_tongue_borders(i) = tmp_S.arnold_tongue_border;
%                 Ntrials(i) = tmp_S.Ntrials;
                score(i) = tmp_S.score;
                score_mean(i) = tmp_S.score_mean;
                score_std(i) = tmp_S.score_std;
            end
        end

        S.arnold_tongue_borders = arnold_tongue_borders;
        S.score = score;
        S.score_mean = score_mean;
        S.score_std = score_std;

        date_string = datestr(clock());
        filename = [output_folder, 'VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=', num2str(S.volume), '_population=', num2str(S.population_average), '_', date_string, '.mat'];
        save(filename, '-struct', 'S');

%         if ~error_occured
%             for i=1:length(input_periods)
%                 filename = [output_folder, 'job_', int2str(n), '.mat'];
%                 display(['deleting file: ', filename]);
% %                 delete(filename);
%             end
%         end

    else
        %% perform a single job

        i = n;

        display(['i=', int2str(i), ' out of ', int2str(length(S.input_periods))]);
        input_period = S.input_periods(i);

        display(['input_period=', num2str(input_period)]);

        lower_amplitude = S.min_input_amplitude;
        upper_amplitude = S.max_input_amplitude;

        upper_amp_scores = simulate_and_compute_entrainment_score_(input_period, upper_amplitude, S.population_average, S);
        upper_amp_within_at = is_within_arnold_tongue__(upper_amp_scores, S);
        lower_amp_scores = simulate_and_compute_entrainment_score_(input_period, lower_amplitude, S.population_average, S);
        lower_amp_within_at = is_within_arnold_tongue__(lower_amp_scores, S);

        mean_upper_amp_score = mean(upper_amp_scores);
        mean_lower_amp_score = mean(lower_amp_scores);

        display(['lower_amplitude=', num2str(lower_amplitude), ', score=', num2str(mean_lower_amp_score), ', within_at=', num2str(lower_amp_within_at)]);
        display(['upper_amplitude=', num2str(upper_amplitude), ', score=', num2str(mean_upper_amp_score), ', within_at=', num2str(upper_amp_within_at)]);

        score = inf;
        score_mean = inf;
        score_std = inf;

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
                scores = simulate_and_compute_all_entrainment_scores_(input_period, middle_amplitude, S.population_average, S);
                middle_amp_within_at = is_within_arnold_tongue__(scores, S);
                mean_score = mean(scores);
                display(['upper_amp=', num2str(upper_amplitude), ', lower_amp=', num2str(lower_amplitude), ', middle_amp=', num2str(middle_amplitude)]);
                display(['  score=', num2str(mean_score), ', within_at=', num2str(middle_amp_within_at)]);
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
            scores = zeros(3, S.Ntrials);
            for j=1:3
                scores(j,:) = simulate_and_compute_entrainment_score_(input_period, input_amplitude, S.population_average, S);
                mean_score = mean(scores(j,:));
                display([' j=', int2str(j), ', score=', num2str(mean_score)]);
            end
            score_mean = mean(scores, 1);
            score_std = std(scores, 1);
            if S.Ntrials == 1
                display(['score_mean=', num2str(score_mean), ', score_std=', num2str(score_std)]);
            end

        end

        display(['arnold tongue border at ', num2str(arnold_tongue_border)]);

        S = struct();
        S.arnold_tongue_border = arnold_tongue_border;
        S.score = score;
        S.scores = scores;
        S.score_mean = score_mean;
        S.score_std = score_std;

        filename = get_job_filename(output_folder, n);
        save(filename, '-struct', 'S');

    end

end

% function filename = get_Ntrials_std_filename(prefix)
%     filename = [prefix, 'Ntrials_std.mat'];
% end

function filename = get_info_filename(folder)
    filename = [folder, 'info.mat'];
end

function filename = get_job_filename(folder, n)
    filename = [folder, 'job_', int2str(n), '.mat'];
end

function within_arnold_tongue = is_within_arnold_tongue__(scores, options)
    if any(isnan(score))
        within_arnold_tongue = true;
    else
        if length(scores) > 1
            ratio_within_arnold_tongue = sum(scores >= options.ENTRAINMENT_THRESHOLD) / length(scores);
            within_arnold_tongue = ratio_within_arnold_tongue >= options.MINIMUM_ENTRAINMENT_RATIO;
        else
            within_arnold_tongue = scores >= options.ENTRAINMENT_THRESHOLD;
        end
    end;
end

% function within_arnold_tongue = is_within_arnold_tongue_(input_period, input_amplitude, population_average, options)
%     if population_average
%         within_arnold_tongue = is_average_within_arnold_tongue(input_period, input_amplitude, options);
%     else
%         within_arnold_tongue = is_within_arnold_tongue(input_period, input_amplitude, options);
%     end
% end

function scores = simulate_and_compute_all_entrainment_scores_(input_period, input_amplitude, population_average, options)
    if population_average
        scores = simulate_average_and_compute_entrainment_scores(input_period, input_amplitude, options);
    else
        scores = simulate_and_compute_all_entrainment_scores(input_period, input_amplitude, options);
    end
end
