addpath('..');
natural_period = 23.8607;

t0 = 0.0;
tf = 200*24;
recordStep = (tf - t0) / 5000.0;
to = (tf - t0) / 10;

input_offset = 1.0;
input_amplitude = 0.3;
input_period = 24.5;
initial_phase = 0;

[T, X] = RunODE(t0, tf, recordStep, ...
    input_offset, input_amplitude, input_period, initial_phase);

%% cutoff transients
offset_time = to;
offset = find(T >= offset_time, 1);
T = T(offset:end);
X = X(offset:end, :);

figure();
plot(T, X);

% [omega, y] = compute_normalized_fft_truncated( ...
%     X(:, 1) - mean(X(:, 1)), recordStep, ...
%     2*pi / 72, 2*pi / 6);
% 
% figure();
% plot(omega / (2*pi), abs(y).^2);


output = X;
min_frequency = 0.0;
max_frequency = 1 / 3;

ENTRAINMENT_THRESHOLD = 0.9;
MAX_HARMONIC_N = 4;
MIN_HARMONICS_POWER_THRESHOLD = 0.0;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.03;
entrainment_ratios = 1:2;

S = struct();
S.natural_period = natural_period;
S.ENTRAINMENT_THRESHOLD = ENTRAINMENT_THRESHOLD;
S.MAX_HARMONIC_N = MAX_HARMONIC_N;
S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
S.entrainment_ratios = entrainment_ratios;


%% substract mean
output = output - repmat(mean(output, 1), [size(output, 1), 1]);

%% Fourier spectrum analysis
[omega, y] = compute_normalized_fft_truncated(output', recordStep, 2*pi*min_frequency, 2*pi*max_frequency);

figure();
plot(omega / (2*pi), abs(y).^2);

score = compute_entrainment_score(omega, y, input_period, S);
score
