natural_period = 2.1013;
ENTRAINMENT_THRESHOLD = 0.9;
MIN_HARMONICS_POWER_THRESHOLD = 0.0;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
% MAX_HARMONIC_N = 15;
MAX_HARMONIC_N = double(intmax());
% entrainment_ratios = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
entrainment_ratios = 1:2;

volume = inf;
omega = volume;

if volume == inf
    Ntrials = 1;
    dt = 0.0001;
    recordStep = 100 * dt;
else
    Ntrials = 100;
    dt = 0.0001;
    recordStep = 100 * dt;
end

t0 = 0;
tf = 1000;
to = (tf - t0) / 5;


input_offset = 1.0;

input_periods = 1.0:0.05:5.0;
input_amplitudes = 0.0:0.05:2.0;

min_frequency = 0.0;
max_frequency = 50.0;

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
S.input_periods = input_periods;
S.input_amplitudes = input_amplitudes;
S.min_frequency = min_frequency;
S.max_frequency = max_frequency;
S.ENTRAINMENT_THRESHOLD = ENTRAINMENT_THRESHOLD;
S.MAX_HARMONIC_N = MAX_HARMONIC_N;
S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
S.entrainment_ratios = entrainment_ratios;

scores = zeros(length(input_periods), length(input_amplitudes));

% for i=1:length(input_periods)
parfor i=1:length(input_periods)
    display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
    input_period = input_periods(i);

    par_scores = zeros(length(input_amplitudes), 1);

    for j=1:length(input_amplitudes)
        display(['i=', int2str(i), ' of ', int2str(length(input_periods)), ', j=', int2str(j), ' of ', int2str(length(input_amplitudes))]);
        input_amplitude = input_amplitudes(j);

        score = simulate_and_compute_entrainment_scores(input_period, input_amplitude, S);
        par_scores(j) = score;
    end

    scores(i, :) = par_scores;

end

S.scores = scores;


date_string = datestr(clock());
filename = ['output/NFkB_TNF_ArnoldTongue_', date_string, '.mat'];
save(filename, '-struct', 'S');
