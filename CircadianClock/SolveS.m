doPlot = true;

t0 = 0;
tf = 20*72;
% um = 1e-6;
um = 1.0;
% x0 = zeros(5, 1);
% x0 = [1.9; 0.8; 0.8; 0.8; 0.8] * um;
x0 = [0.1; 0.25; 0.25; 0.25; 0.25] * um;

dt = 0.001;
recordStep = 100 * dt;
Nsteps = ceil(tf / dt);
Ntrials = 50;

%% Results
% Amplitude 0.2 period 36
% - for volume inf natural mode has a power of 5e7, input mode has a power
% of 3e7
% - effect seen for volume 1e-21, reduction of natural mode from 2e7 to
%   0.5e7, input mode remains at 3e7
% - effect seen for volume 5e-22, reduction of natural mode from ~2e7 to
%   ~0.2e7, input mode remains at 3e7
% - effect seen for volume 1e-22, reduction of natural mode from ~1e7 to
%   <0.1e7, input mode remains at 3.5e7
%
% Amplitude 0.15 period 36
% - for volume inf natural mode has a power of 8e7, input mode has a power
% of 2e7
% - effect seen for volume 5e-22, reduction of natural mode from 3e7 to
%   ~0.4e7, input mode remains at ~2e7



% natural_period = 23.6574;
% input_amplitude = 0.0;
input_amplitude = 0.2;
input_offset = 1.0;
input_period = 36.0;
input_frequency = 1 ./ input_period;
input_function = @(t, x) input_offset + input_amplitude * sin(2 * pi * t ./ input_period);

% volume = inf;
volume = 1e-20;
volume_str = ['volume=', num2str(volume), '> '];

% seed = 858526701;
% [T, X, omega] = SolveS_Java(x0, tf, dt, volume, ...
%                             input_offset, input_amplitude, input_frequency, Ntrials, ...
%                             false, recordStep, seed);
[T, X, omega] = SolveS_Java(x0, tf, dt, volume, ...
                            input_offset, input_amplitude, input_frequency, Ntrials, ...
                            false, recordStep);
P = X;
X = squeeze(X(1, :, :));

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
[omega, y] = compute_fft(cutoff_P_t, recordStep);
freq = omega ./ (2 * pi);
min_freq = (1 / 24) / 2;
max_freq = (1 / 24) * 2;
i1 = find(freq >= min_freq, 1);
i2 = find(freq <= max_freq, 1, 'last');
figure();
plot(freq(i1:i2), abs(y(i1:i2)).^2);


Nn = sum(P(:, :, 2:5), 3);

omega = [];
y = [];
for i=1:Ntrials
    [omega1, y1]= compute_fft(Nn(i,:), recordStep);
    min_omega = 2*pi * (1 / 24) / 2;
    max_omega = 2*pi * (1 / 24) * 2;
    i1 = find(omega1 < min_omega, 1, 'last');
    if length(i1) == 0
        i1 = 1;
    end
    i2 = find(omega1 > max_omega, 1, 'first');
    if length(i2) == 0
        i2 = length(omega1);
    end
    omega1 = omega1(i1:i2);
    y1 = y1(i1:i2);
    omega = [omega; omega1];
    y = [y; y1];
end

mean_omega = mean(omega, 1);
mean_freq = mean(omega ./ (2 * pi), 1);
% mean_y = mean(y, 1);
abs_mean_omega = abs(mean(y, 1)) .^ 2;
mean_abs_omega = mean(abs(y) .^ 2);

dfreq = 0.05;
disp('average around frequency 0.33');
freq = 0.33;
freq1 = freq;
mask = (mean_freq >= freq-dfreq) & (mean_freq <= freq+dfreq);
q1 = sum(abs_mean_omega(mask));
disp([' complex average:', num2str(q1)]);
if Ntrials > 1
    w1 = sum(mean_abs_omega(mask));
    disp([' absolute average:', num2str(w1)]);
end
disp('average around frequency 0.50');
freq = 0.50;
freq2 = freq;
mask = (mean_freq >= freq-dfreq) & (mean_freq <= freq+dfreq);
q2 = sum(abs_mean_omega(mask));
disp([' complex average:', num2str(q2)]);
if Ntrials > 1
    w2 = sum(mean_abs_omega(mask));
    disp([' absolute average:', num2str(w2)]);
end
disp(['ratios ', num2str(freq2), ' to ', num2str(freq1)]);
disp([' complex: ', num2str(q2 / q1)]);
if Ntrials > 1
    disp([' absolute: ', num2str(w2 / w1)]);
end

if doPlot
    figure;
    subplot(2, 1, 1);
    plot(mean_freq, abs_mean_omega);
    title(['complex average spectrum for volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
    xlabel('frequency f');
    ylabel('|y|^{2}');
%     ylim([0, 1e6]);
    subplot(2, 1, 2);
    plot(mean_freq, mean_abs_omega);
    title(['absolute average spectrum for volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
    xlabel('frequency f');
    ylabel('|y|^{2}');
end
