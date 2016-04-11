addpath('..');
natural_period = 24.5761;

% volume = 1.7e-15;
% volume = 1 / 6.022e23;
volume = 9.9635e-22;
volume = volume * 1;
% omega = inf;

Ntrials = 1;

t0 = 0;
dt = 0.01;
tf = 10*24;
recordStep = (tf - t0)/1000;

[T, X, omega, model] = Run(Ntrials, t0, tf, dt, recordStep, volume, ...
    1.0, 0.0, 1.1 * natural_period, 0.0);
[omega, y] = compute_normalized_fft_truncated( ...
    X(:, 1) - mean(X(:, 1)), recordStep, ...
    2*pi / 72, 2*pi / 6);

prop = model.getSSAModel().computePropensities(0, X(end,:))

figure();
plot(T, X(:, 1));

figure();
plot(omega / (2*pi), abs(y).^2);
