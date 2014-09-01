% brutus
% POPULATION_AVERAGE=false
% POPULATION_AVERAGE=true
% bsub -n 1 -R "rusage[mem=1536]" -W 12:00 -J "job[1-525]" -o logs2/CircadianClock_ArnoldTongue_JobArray_%I.out bash CircadianClock_ArnoldTongue_JobArray.sh "\$LSB_JOBINDEX" output2/ $POPULATION_AVERAGE
% bsub -n 1 -R "rusage[mem=1536]" -W 8:00 -J "job_comb" -o logs2/CircadianClock_ArnoldTongue_JobArray_Combine.out bash CircadianClock_ArnoldTongue_JobArray.sh 0 output2/ $POPULATION_AVERAGE
%
% INDEX=71; bsub -n 1 -R "rusage[mem=2048]" -W 16:00 -o logs2/CircadianClock_ArnoldTongue_JobArray_$INDEX.out bash CircadianClock_ArnoldTongue_JobArray.sh $INDEX output2/ $POPULATION_AVERAGE
%

function CircadianClock_ArnoldTongue_JobArray(n, filename_prefix, population_average)

    if nargin < 3
        population_average = false;
    end

    ENTRAINMENT_THRESHOLD = 0.9;
    MAX_HARMONIC_N = 4;
    MIN_HARMONICS_POWER_THRESHOLD = 1.0;
    FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
%     STD_ESTIMATION_SIZE = 3;
    natural_period = 23.7473;

%     volume = inf;
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
    tf = 10000;
    to = (tf - t0) / 5;


    input_offset = 1.0;

    input_periods = 12:1.0:36;
    input_amplitudes = 0.0:0.05:1.0;

    min_frequency = 0.005;
    max_frequency = 0.5;

    [input_period_mesh, input_amplitude_mesh] = meshgrid(input_periods, input_amplitudes);


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
    S.natural_period = natural_period;
    S.min_frequency = min_frequency;
    S.max_frequency = max_frequency;
    S.input_offset = input_offset;
    S.input_periods = input_periods;
    S.input_amplitudes = input_amplitudes;
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

        scores = zeros(length(input_periods), length(input_amplitudes));

        for i=1:length(input_periods)
            display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
            for j=1:length(input_amplitudes)
                display(['j=', int2str(j), ' out of ', int2str(length(input_amplitudes))]);
                id = (i- 1) * length(input_amplitudes) + j;
                filename = get_filename(filename_prefix, id);
                error_occured = false;
                try
                    tmp_S = load(filename);
                catch
                    error_occured = true;
                    display(['error occured: skipping i=', int2str(i), ' j=', int2str(j)]);
                end
                if ~error_occured
                    scores(i, j) = tmp_S.score;
                end
            end
        end

        S.natural_period = natural_period;
        S.input_periods = input_periods;
        S.scores = scores;

        date_string = datestr(clock());
        filename = [output_folder, 'CircadianClock_ArnoldTongue_JobArray_volume=', num2str(volume), '_', date_string, '.mat'];
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

        display(['i=', int2str(i), ' out of ', int2str(length(input_periods) * length(input_amplitudes))]);
        input_period = input_period_mesh(i);
        input_amplitude = input_amplitude_mesh(i);

        score = simulate_and_compute_entrainment_score_(input_period, input_amplitude, population_average, S);

        display(['i=', int2str(i), ' score=', num2str(score)]);

        S = struct();
        S.score = score;

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

function score = simulate_and_compute_entrainment_score_(input_period, input_amplitude, population_average, options)
    if population_average
        score = simulate_average_and_compute_entrainment_scores(input_period, input_amplitude, options);
    else
        score = simulate_and_compute_entrainment_scores(input_period, input_amplitude, options);
    end
end
