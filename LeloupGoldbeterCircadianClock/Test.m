addpath('..');
natural_period = 24.5761;

Ntrials = 8;

omega = 600;
t0 = 0.0
tf = 20 * 24;
recordStep = tf/1000;

[T, X, model] = Run(Ntrials, t0, tf, recordStep, omega, ...
    1.0, 0.0, 1.1*natural_period, 0.0);

figure();
plot(T, X);

[omega, y] = compute_normalized_fft_truncated( ...
    X(1, :) - mean(X(1, :)), recordStep, ...
    2*pi / 72, 2*pi / 6);
figure();
plot(omega / (2*pi), abs(y).^2);
