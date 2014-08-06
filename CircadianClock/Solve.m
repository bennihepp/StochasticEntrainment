addpath('..');

model = @CircadianClock;
t0 = 0;
tf = 20*72;
% um = 1e-6;
um = 1.0;
x0 = zeros(5, 1);
% x0 = [1.9; 0.8; 0.8; 0.8; 0.8] * um;
% x0 = [0.1; 0.25; 0.25; 0.25; 0.25] * um;

% natural_period = 23.6574;
% input_amplitude = 0.0;
input_amplitude = 0.0;
input_offset = 1.0;
input_period = 28.0;
input_function = @(t, x) input_offset + input_amplitude * sin(2 * pi * t ./ input_period);

% tspan = [0, tf];
dt = 0.1;
tspan = t0:dt:tf;
[T, X] = ode45(model, tspan, x0, [], input_function);

figure();
M = X(:, 1);
P_0 = X(:, 2);
P_t = sum(X(:, 2:5), 2);
inp = input_function(T);
plot(T, M, 'b');
hold on;
plot(T, P_0, 'g');
plot(T, P_t, 'c');
plot(T, inp, 'r');
legend('M', 'P_0', 'P_t', 'input');

cutoff_P_t = P_t(floor(length(P_t) / 2):end);
[omega, y] = compute_fft(cutoff_P_t, dt);
freq = omega ./ (2 * pi);
min_freq = (1 / 24) / 2;
max_freq = (1 / 24) * 2;
i1 = find(freq >= min_freq, 1);
i2 = find(freq <= max_freq, 1, 'last');
figure();
plot(freq(i1:i2), abs(y(i1:i2)).^2);
