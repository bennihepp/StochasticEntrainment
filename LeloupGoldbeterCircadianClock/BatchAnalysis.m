%% plot trajectories
figure();
plot(T, output(:, 1));
title(['y(1) first trace: Ntrials=', int2str(Ntrials), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('time t');
ylabel('state y(1)');

figure();
plot(T, mean(output, 2));
title(['y(1) average trace: Ntrials=', int2str(Ntrials), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('time t');
ylabel('state y(1)');


%% cutoff transients
% offset_time = (tf - t0) / 5;
% offset_time = max(offset_time, 1000);
offset_time = to;
offset = find(T >= offset_time, 1);
T = T(offset:end);
output = output(offset:end, :);


%% substract mean
output = output - repmat(mean(output, 1), [size(output, 1), 1]);

%% compute spectras
addpath('../');

% min_frequency = 0.005;
% max_frequency = 0.5;
min_frequency = 0.0;
max_frequency = 1 / 6;

omega = [];
y = [];
for i=Ntrials:-1:1
    display([int2str(i), ' out of ', int2str(Ntrials)]);
    [omega1, y1] = compute_normalized_fft_truncated(output(:,i)', recordStep, 2*pi*min_frequency, 2*pi*max_frequency);
    omega = [omega; omega1];
    y = [y; y1];
end
mean_y = mean(y, 1);
mean_omega = mean(omega, 1);

%% plot spectras
figure();
plot(mean_omega ./ (2 * pi), mean(abs(y(1,:)), 1) .^ 2);
title(['y(1) first trace fft: Ntrials=', int2str(Ntrials), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('frequency f');
ylabel('power |y|^2');

figure();
plot(mean_omega ./ (2 * pi), abs(mean_y) .^ 2);
title(['y(1) complex average fft: Ntrials=', int2str(Ntrials), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('frequency f');
ylabel('power |y|^2');

figure();
plot(mean_omega ./ (2 * pi), mean(abs(y), 1) .^ 2);
title(['y(1) absolute average fft: Ntrials=', int2str(Ntrials), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('frequency f');
ylabel('power |y|^2');


%% plot phase distribution of natural mode and input mode
NUM_OF_BINS = 100;
bins = linspace(-pi, pi, NUM_OF_BINS);
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ natural_period));
figure();
hist(angle(y(:, ind)), bins);
title(['phase distribution of natural mode for input offset=', num2str(input_offset), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');
% saveas(['phase_distribution_natural_mode_volume=', num2str(volume), '_input_period=', num2str(input_period), '_input_amplitude=', num2str(input_amplitude), '_Ntrials=', int2str(Ntrials), '.fig']);

[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
figure();
hist(angle(y(:, ind)), bins);
title(['phase distribution of input mode for input offset=', num2str(input_offset), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');
% saveas(['phase_distribution_input_mode_volume=', num2str(volume), '_input_period=', num2str(input_period), '_input_amplitude=', num2str(input_amplitude), '_Ntrials=', int2str(Ntrials), '.fig']);
