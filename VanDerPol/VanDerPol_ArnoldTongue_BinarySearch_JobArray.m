% brutus
% POPULATION_AVERAGE=false
% POPULATION_AVERAGE=true
% bsub -n 1 -R "rusage[mem=1536]" -W 12:00 -J "job[1-81]" -o logs2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_%I.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh "\$LSB_JOBINDEX" output2/ $POPULATION_AVERAGE
% bsub -n 1 -R "rusage[mem=1536]" -W 8:00 -J "job_comb" -o logs2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_Combine.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh 0 output2/
%
% INDEX=71; bsub -n 1 -R "rusage[mem=2048]" -W 16:00 -o logs2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_$INDEX.out bash VanDerPol_ArnoldTongue_BinarySearch_JobArray.sh $INDEX output2/ $POPULATION_AVERAGE
%

function VanDerPol_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, population_average)

    if nargin < 3
        population_average = false;
    end

    ENTRAINMENT_THRESHOLD = 0.9;
    MAX_HARMONIC_N = 4;
    MIN_HARMONICS_POWER_THRESHOLD = 1.0;
    FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
    natural_period = 1/0.1065;

    volume = 5e3;

    if volume == inf
        Ntrials = 1;
        dt = 1e-1;
    else
        Ntrials = 100;
        dt = 1e-1;
    end

    disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);

    t0 = 0;
    tf = 10000;
    to = (tf - t0) / 5;


    min_input_amplitude = 0.0;
    max_input_amplitude = 2.0;
    input_amplitude_tolerance = 2e-2;
    input_periods = 4:0.2:20.0;

    min_frequency = 0.01;
    max_frequency = 1.0;


    S = struct();
    S.volume = volume;
    S.Ntrials = Ntrials;
    S.t0 = t0;
    S.tf = tf;
    S.to = to;
    S.dt = dt;
    S.natural_period = natural_period;
    S.min_frequency = min_frequency;
    S.max_frequency = max_frequency;
    S.ENTRAINMENT_THRESHOLD = ENTRAINMENT_THRESHOLD;
    S.MAX_HARMONIC_N = MAX_HARMONIC_N;
    S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
    S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;


    if n == 0

        arnold_tongue_borders = zeros(length(input_periods), 1);
        score_variances = zeros(length(input_periods), 1);

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
                score_variances(i) = tmp_S.score_variance;
            end
        end

        S.natural_period = natural_period;
        S.input_periods = input_periods;
        S.min_input_amplitude = min_input_amplitude;
        S.max_input_amplitude = max_input_amplitude;
        S.input_amplitude_tolerance = input_amplitude_tolerance;
        S.arnold_tongue_borders = arnold_tongue_borders;
        S.score_variances = score_variances;

        date_string = datestr(clock());
        filename = ['VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=', num2str(volume), '_', date_string, '.mat'];
        save(filename, '-struct', 'S');

%         if ~error_occured
%             for i=1:length(input_periods)
%                 filename = [filename_prefix, 'job_', int2str(n), '.mat'];
%                 display(['deleting file: ', filename]);
% %                 delete(filename);
%             end
%         end

    else

        i = n;

        display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
        input_period = input_periods(i);

        lower_amplitude = min_input_amplitude;
        upper_amplitude = max_input_amplitude;

        upper_amp_within_at = is_within_arnold_tongue_(input_period, upper_amplitude, population_average, S);
        lower_amp_within_at = is_within_arnold_tongue_(input_period, lower_amplitude, population_average, S);

        score_variance = inf;

        if lower_amp_within_at
            arnold_tongue_border = lower_amplitude;
        elseif ~upper_amp_within_at
            arnold_tongue_border = inf;
        else

            while (upper_amplitude - lower_amplitude) >= input_amplitude_tolerance
                middle_amplitude = (lower_amplitude + upper_amplitude) / 2.0;
                display(['i=', int2str(i), ' of ', int2str(length(input_periods)), ', trying input=', num2str(middle_amplitude)]);
                middle_amp_within_at = is_within_arnold_tongue_(input_period, middle_amplitude, population_average, S);
                if middle_amp_within_at
                    upper_amplitude = middle_amplitude;
                else
                    lower_amplitude = middle_amplitude;
                end
            end

            if middle_amp_within_at
                arnold_tongue_border = middle_amplitude;
            else
                arnold_tongue_border = upper_amplitude;
            end

            scores = zeros(3, 1);
            input_amplitude = middle_amplitude;
            for j=1:3
                scores(j) = simulate_and_compute_entrainment_scores(input_period, input_amplitude, S);
            end
            score_variance = std(scores);

        end

        display(['i=', int2str(i), ' arnold tongue border at ', num2str(arnold_tongue_border)]);

        S = struct();
        S.arnold_tongue_border = arnold_tongue_border;
        S.score_variance = score_variance;

        filename = get_filename(filename_prefix, n);
        save(filename, '-struct', 'S');

    end

end

function filename = get_filename(prefix, n)
    filename = [prefix, 'job_', int2str(n), '.mat'];
end

function within_arnold_tongue = is_within_arnold_tongue_(input_period, input_amplitude, population_average, options)
    if population_average
        within_arnold_tongue = is_average_within_arnold_tongue(input_period, input_amplitude, options);
    else
        within_arnold_tongue = is_within_arnold_tongue(input_period, input_amplitude, options);
    end
end
