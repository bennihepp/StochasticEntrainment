%% some functions in the parent folder are used
addpath('../');
addpath('../plotting/');

natural_period = 23.8607;

t0 = 0.0;
tf = 400*24;
recordStep = (tf - t0) / 5000.0;
to = (tf - t0) / 10;

input_offset = 1.0;
initial_phase = 0;

% input_period = 23.5;
% input_amplitude = 0.3;

% input_period = 22.2;
% input_amplitude = 0.25;

input_period = 24.0;
input_amplitude = 0.2;

input_period = 26.0;
input_amplitude = 0.2;

to = 600;
recordStep = 0.12;
tf = 1200;

min_frequency = 0.0;
max_frequency = inf;

[T, X] = RunODE(t0, tf, recordStep, ...
    input_offset, input_amplitude, input_period, initial_phase);

%% cutoff transients
offset_time = to;
offset = find(T >= offset_time, 1);
T = T(offset:end);
X = X(offset:end, :);

return;

%% plot traces after transients

output = X;
w = find(T > T(end) - 400, 1);
q = find(T > T(w) - 150, 1);
TT = T(q:w);
TT = TT - TT(1);
trunc = output(q:w, :);

width = 10;
height = 3;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
hold on;
deterministic_color = [0.9, 0.0, 0.0];
plot(TT, mean(trunc, 2), '-', 'Color', deterministic_color, 'LineWidth', 1.0);
xlabel('time t');
ylabel('state y');
hold off;
% save_plot([export_eps_prefix(), 'leloup_goldbeter_circadian_deterministic_trace'], h, width, height);


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
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
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
