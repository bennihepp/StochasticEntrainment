PERIODICITY_THRESHOLD = 0.01;

volume = inf;
% volume = 1e3;
% volume = 1e2;
% volume = 5e1;
% volume = 1e1;
% volume = 1e0;
% volume = 1e-1;
disp(volume);

omega = volume;

if volume == inf
    Ntrials = 1;
    dt = 1e-1;
else
    Ntrials = 50;
    dt = 1e-1;
end

t0 = 0
tf = 10000;

% % input_amplitude = 1.0;
% % input_amplitude = 0.8;
% input_amplitude = 0.0;
% input_period = 5;

input_period = 19;
input_amplitude = 0.5;

additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
multiplicative_forcing_func = @(t, x) 0;


[T, output] = VanDerPol_Run(Ntrials, t0, tf, dt, omega, additive_forcing_func, multiplicative_forcing_func);


figure();
plot(T, output(:, 1));
title(['y(1) first trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('time t');
ylabel('state y(1)');

figure();
plot(T, mean(output, 2));
title(['y(1) average trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('time t');
ylabel('state y(1)');


offset_time = (tf - t0) / 5;
offset_time = min(offset_time, 1000);
offset = find(T >= offset_time, 1);
T = T(offset:end);
output = output(offset:end, :);


addpath('../');

omega = [];
y = [];
for i=Ntrials:-1:1
    [omega1, y1]= compute_fft(output(:,i)', dt);
    min_omega = 2 * pi * 0.01;
    max_omega = 2 * pi * 1.0;
    i1 = find(omega1 < min_omega, 1, 'last');
    i2 = find(omega1 > max_omega, 1, 'first');
    omega1 = omega1(i1:i2);
    y1 = y1(i1:i2);
    omega = [omega; omega1];
    y = [y; y1];
end
mean_y = mean(y, 1);
mean_omega = mean(omega, 1);

figure();
plot(mean_omega ./ (2 * pi), abs(mean_y) .^ 2);
title(['y(1) complex average fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('frequency f');
ylabel('power |y|^2');

figure();
plot(mean_omega ./ (2 * pi), mean(abs(y), 1) .^ 2);
title(['y(1) absolute average fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('frequency f');
ylabel('power |y|^2');

% filename = ['output/simulation_Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' offset=', num2str(TNF_offset), ' amplitude=', num2str(TNF_amplitude), ' period=', num2str(TNF_period), '.mat'];
% save(filename);

corr = xcorr(output - mean(output), 'unbiased');
figure();
plot(corr);
[pks, locs] = findpeaks(corr);
peak_distances = locs(2:end) - locs(1:end-1);
mean_peak_distance = mean(peak_distances);
std_peak_distance = std(peak_distances);

output_period = mean_peak_distance * dt;
output_periodic = false;
if std_peak_distance / mean_peak_distance < PERIODICITY_THRESHOLD
    output_periodic = true;
end
