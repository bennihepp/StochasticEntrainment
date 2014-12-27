addpath('..');
natural_period = 23.8607;

t0 = 0.0;
tf = 100*24;
recordStep = (tf - t0) / 1000.0;

input_offset = 1.0;
input_amplitude = 0.0;
input_period = 23.8607;
initial_phase = 0;

[T, X] = RunODE(t0, tf, recordStep, ...
    input_offset, input_amplitude, input_period, initial_phase);
[omega, y] = compute_normalized_fft_truncated( ...
    X(:, 1) - mean(X(:, 1)), recordStep, ...
    2*pi / 72, 2*pi / 6);

figure();
plot(T, X(:, 1));

figure();
plot(omega / (2*pi), abs(y).^2);
