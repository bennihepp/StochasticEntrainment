ENTRAINMENT_THRESHOLD = 0.9;
MAX_HARMONIC_N = 4;
MIN_HARMONICS_POWER_THRESHOLD = 1.0;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
natural_period = 23.7473;

volume = inf;
omega = volume;

if volume == inf
    Ntrials = 1;
    dt = 0.002;
    recordStep = 100 * dt;
else
    Ntrials = 100;
    dt = 0.002;
    recordStep = 100 * dt;
end

t0 = 0;
tf = 100*72;
to = (tf - t0) / 5;


input_periods = 12:1.0:36;
input_amplitudes = 0.0:0.05:1.0;

min_frequency = 0.005;
max_frequency = 0.5;


S = struct();
S.volume = volume;
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
S.ENTRAINMENT_THRESHOLD = ENTRAINMENT_THRESHOLD;
S.MAX_HARMONIC_N = MAX_HARMONIC_N;
S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;


for i=1:length(input_periods)
% parfor i=1:length(input_periods)
    display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
    input_period = input_periods(i);

    lower_amplitude = min_input_amplitude;
    upper_amplitude = max_input_amplitude;

    upper_amp_within_at = is_within_arnold_tongue(input_period, upper_amplitude, S);
    lower_amp_within_at = is_within_arnold_tongue(input_period, lower_amplitude, S);

    if lower_amp_within_at
        arnold_tongue_borders(i) = lower_amplitude;
    elseif ~upper_amp_within_at
        arnold_tongue_borders(i) = inf;
    else

        while (upper_amplitude - lower_amplitude) >= input_amplitude_tolerance
            middle_amplitude = (lower_amplitude + upper_amplitude) / 2.0;
            display(['i=', int2str(i), ' of ', int2str(length(input_periods)), ', trying input=', num2str(middle_amplitude)]);
            middle_amp_within_at = is_within_arnold_tongue(input_period, middle_amplitude, S);
            if middle_amp_within_at
                upper_amplitude = middle_amplitude;
            else
                lower_amplitude = middle_amplitude;
            end
        end

        if middle_amp_within_at
            arnold_tongue_borders(i) = middle_amplitude;
        else
            arnold_tongue_borders(i) = upper_amplitude;
        end

    end

    if isinf(arnold_tongue_borders(i))
        scores(i) = nan;
    else
        scores(i) = simulate_and_compute_entrainment_scores(input_period, arnold_tongue_borders(i), S);
    end

    display(['i=', int2str(i), ' arnold tongue border at ', num2str(arnold_tongue_borders(i))]);

end

S.natural_period = natural_period;
S.input_periods = input_periods;
S.min_input_amplitude = min_input_amplitude;
S.max_input_amplitude = max_input_amplitude;
S.input_amplitude_tolerance = input_amplitude_tolerance;
S.arnold_tongue_borders = arnold_tongue_borders;
S.scores = scores;


date_string = datestr(clock());
filename = ['CircadianClock_ArnoldTongue_BinarySearch_', date_string, '.mat'];
save(filename, '-struct', 'S');
