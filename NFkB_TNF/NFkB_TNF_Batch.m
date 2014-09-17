natural_period = 2.1013;
PERIOD_DEVIATION_THRESHOLD = 0.01 * natural_period;
PERIODICITY_THRESHOLD = 0.05;
PERIOD_MULTIPLE_THRESHOLD = 0.01;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.05;
MIN_HARMONICS_POWER_THRESHOLD = 0.0;
% MAX_HARMONIC_N = 15;
MAX_HARMONIC_N = double(intmax());
% entrainment_ratios = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
entrainment_ratios = 1:2;


%% results for input_period = 2.28, input_amplitude = 0.016

% volume = 1e-10; % with Ntrials=500, population=0.92258, individual=0.10358
% volume = 5e-11; % with Ntrials=500, population=0.92245, individual=0.21997
% volume = 2e-11; % with Ntrials=500, population=0.9321, individual=0.35831
% volume = 1e-11; % with Ntrials=500, population=0.93519, individual=0.32637
% volume = 5e-12; % with Ntrials=500, population=0.93405, individual=0.27137
% volume = 1e-12; % with Ntrials=500, population=0.90784, individual=0.14164

% volume = 5e-11; % with Ntrials=100, population=0.88063, individual=0.20546
% volume = 2e-11; % with Ntrials=100, population=0.88031, individual=0.3355
% volume = 1e-11; % with Ntrials=100, population=0.8837, individual=0.32711
% volume = 5e-12; % with Ntrials=100, population=0.879, individual=0.27711
% volume = 1e-12; % with Ntrials=100, population=0.80449, individual=0.1464
% volume = 5e-13; % with Ntrials=100, population=0.71225, individual=0.10105
% volume = 1e-13; % with Ntrials=100, population=0.18849, individual=0.053604


%% parameters

% volume = inf;
% volume = 1e-10;

% volume = 5e-12;
% volume = 1e-11;
volume = 2e-11;

if volume == inf
    Ntrials = 1;
else
    Ntrials = 100;
    Ntrials = 50;
end

dt = 0.0001;
recordStep = 100 * dt;

disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);

%% complex average entrainment score with Ntrials=100, input_period=2.6, input_amplitude=0.055
% volume=1e-10:  -
% volume=5e-11:  0.80324
% volume=2e-11:  0.82477
% volume=1e-11:  0.80033
% volume=5e-12:  0.71624
% volume=2e-12:  0.74775
% volume=1e-12:  0.76091
% volume=5e-13:  0.87799
% volume=2e-13:  0.86608
% volume=1e-13:  0.80698
% volume=5e-14:  0.63785

%% parameters
t0 = 0;
tf = 1000;
to = (tf - t0) / 5;

input_offset = 1.0;

% input_period = 36.0;
% input_amplitude = 0.5;

input_period = 2.1013;
input_amplitude = 0.2;

input_period = 2.28;
input_amplitude = 0.1;

input_period = 2.28;
input_amplitude = 0.016;
input_amplitude = 0.05;

input_period = 2.4;
input_amplitude = 0.025;

input_period = 2.24;
% input_period = 2.26;
input_amplitude = 0.05;

% min_frequency = 0.001;
min_frequency = 0.0;
% max_frequency = inf;
% max_frequency = 50.0;
max_frequency = 10.0;


%% create options structure
S = struct();
S.min_frequency = min_frequency;
S.max_frequency = max_frequency;
S.natural_period = natural_period;
S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
S.MAX_HARMONIC_N = MAX_HARMONIC_N;
S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
S.entrainment_ratios = entrainment_ratios;


%% simulate
tic;
printProgress = true;
[T, original_output] = NFkB_TNF_Run(Ntrials, t0, tf, dt, recordStep, volume, input_offset, input_amplitude, input_period, printProgress);
toc
output = original_output;

% filename = ['output/simulation_Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), '_offset=', num2str(input_offset), ',_amplitude=', num2str(input_amplitude), ',_period=', num2str(input_period), '.mat'];
% save(filename);


% %% plot trajectories
% figure();
% plot(T, output(:, 1));
% hold on;
% plot([to, to], ylim(), '-r');
% hold off;
% title(['y(1) first trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
% xlabel('time t');
% ylabel('state y(1)');
% 
% figure();
% plot(T, mean(output, 2));
% hold on;
% plot([to, to], ylim(), '-r');
% hold off;
% title(['y(1) average trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
% xlabel('time t');
% ylabel('state y(1)');


%% cutoff transients
% offset_time = (tf - t0) / 5;
% offset_time = max(offset_time, 1000);
offset_time = to;
offset = find(T >= offset_time, 1);
T = T(offset:end);
output = output(offset:end, :);


%% plot traces after transients

w = find(T > T(end) - 50, 1);
q = find(T > T(w) - 10, 1);
TT = T(q:w);
TT = TT - TT(1);
trunc = output(q:w, :);

figure;
%     plot(T, output(:, i), 'Color', cmap(i, :), 'LineWidth', 2.0);
plot(TT, trunc(:, 1), 'LineWidth', 2.0);
title(['y(1) single trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
xlabel('time t');
ylabel('state y');

figure;
hold on;
cmap = colormap('Lines');
for i=1:3
%     plot(T, output(:, i), 'Color', cmap(i, :), 'LineWidth', 2.0);
    plot(TT, trunc(:, i), 'Color', cmap(i, :), 'LineWidth', 2.0);
end
hold off;
title(['y(1) single traces: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
xlabel('time t');
ylabel('state y');

figure();
% plot(T, mean(output, 2), 'LineWidth', 2.0);
plot(TT, mean(trunc, 2), 'LineWidth', 2.0);
title(['y(1) average trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('time t');
ylabel('state y(1)');


%% substract mean
output = output - repmat(mean(output, 1), [size(output, 1), 1]);

%% plot corrected output
figure();
plot(T, mean(output, 2));
hold on;
plot([to, to], ylim(), '-r');
hold off;
title(['y(1) corrected average trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('time t');
ylabel('state y(1)');


%% compute spectras
addpath('../');

omega = [];
y = [];
for i=Ntrials:-1:1
    [omega1, y1] = compute_normalized_fft_truncated(output(:,i)', recordStep, 2*pi*min_frequency, 2*pi*max_frequency);
    omega = [omega; omega1];
    y = [y; y1];
end
mean_y = mean(y, 1);
mean_omega = mean(omega, 1);

%% plot spectras
figure();
plot_spectrum_and_mark_harmonics(mean_omega, y(1, :), 2*pi/input_period, 'red', S, 1);
% figure();
% plot(mean_omega ./ (2 * pi), mean(abs(y(1,:)), 1) .^ 2);
title(['y(1) first trace fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);

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


%% plot phase distribution of natural mode and input mode
NUM_OF_BINS = 100;
bins = linspace(-pi, pi, NUM_OF_BINS);
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ natural_period));
figure();
hist(angle(y(:, ind)), bins);
title(['phase distribution of natural mode for volume=', num2str(volume), ' input offset=', num2str(input_offset), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');
% saveas(['phase_distribution_natural_mode_volume=', num2str(volume), '_input_period=', num2str(input_period), '_input_amplitude=', num2str(input_amplitude), '_Ntrials=', int2str(Ntrials), '.fig']);

[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
figure();
hist(angle(y(:, ind)), bins);
title(['phase distribution of input mode for volume=', num2str(volume), ' input offset=', num2str(input_offset), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');
% saveas(['phase_distribution_input_mode_volume=', num2str(volume), '_input_period=', num2str(input_period), '_input_amplitude=', num2str(input_amplitude), '_Ntrials=', int2str(Ntrials), '.fig']);


%% compute autocorrelation if only one trajectory is simulated
if Ntrials == 1
    corr = xcorr(output - mean(output, 1), 'unbiased');
    figure();
    plot(corr);
    title('autocorrelation');
    [pks, locs] = findpeaks(corr);
    peak_distances = locs(2:end) - locs(1:end-1);
    mean_peak_distance = mean(peak_distances);
    std_peak_distance = std(peak_distances);

    mean_period = mean_peak_distance * recordStep;
    factor = mean_period / input_period;
    if mean_period < input_period
        factor = input_period / mean_period;
    end
    output_periodic = false;
    if abs(factor - round(factor)) < PERIOD_MULTIPLE_THRESHOLD
    % if abs(mean_period - input_period) < PERIOD_DEVIATION_THRESHOLD
        output_periodic = std_peak_distance / mean_peak_distance < PERIODICITY_THRESHOLD;
    end
    output_periodic

%     [q,w] = compute_fft(corr, recordStep);
%     figure();
%     plot(q/2/pi, abs(w).^2);
end


%% compute entrainment
% om_natural = 2 * pi / natural_period;
% om_input = 2 * pi / input_period;
% dom = FREQUENCY_NEIGHBOURHOOD_FACTOR * om_natural;
% 
% power_total = sum(abs(mean_y).^2);
% power_input = compute_spectrum_power(mean_omega, mean_y, om_input, dom);
% power_input_harmonics = 0;
% for n=2:MAX_HARMONIC_N
%     power_input_harmonics = power_input_harmonics + compute_spectrum_power(mean_omega, mean_y, om_input * n, dom);
% end
% if power_input >= 0.1 * power_input_harmonics
%     power_input = power_input + power_input_harmonics;
% end
% power_input / power_total


%% compute entrainment scores
% Omega = mean_omega;
% Q = zeros(Ntrials, 1);
% W = zeros(Ntrials, 1);
% om_natural = 2 * pi / natural_period;
% om_input = 2 * pi / input_period;
% dom = FREQUENCY_NEIGHBOURHOOD_FACTOR * om_natural;
% if abs(om_natural - om_input) < 2 * dom
%     Q(:) = inf;
%     W(:) = inf;
% else
%     for n=1:Ntrials
%         power_total = sum(abs(y(n,:)).^2);
%         power_natural = compute_spectrum_power(Omega, y(n,:), om_natural, dom);
% %         for n=2:MAX_HARMONIC_N
% %             power_natural = power_natural + compute_spectrum_power(Omega, y(n,:), om_natural * n, dom);
% %         end
%         power_input = compute_spectrum_power(Omega, y(n,:), om_input, dom);
%         power_input_harmonics = 0;
%         for m=2:MAX_HARMONIC_N
%             power_input_harmonics = power_input_harmonics + compute_spectrum_power(Omega, y(n,:), om_input * m, dom);
%         end
%         if power_input >= MIN_HARMONICS_POWER_THRESHOLD * power_input_harmonics
%             power_input = power_input + power_input_harmonics;
%         end
% 
%         Q(n) = power_input / power_natural;
%         W(n) = power_input / power_total;
% 
%     end
% 
% end

% if abs(om_natural - om_input) < 2 * dom
%     Q_mean = inf;
%     W_mean = inf;
% else
%     power_total = sum(abs(mean_y).^2);
%     power_natural = compute_spectrum_power(Omega, mean_y, om_natural, dom);
% %     for n=2:MAX_HARMONIC_N
% %         power_natural = power_natural + compute_spectrum_power(Omega, mean_y, om_natural * n, dom);
% %     end
%     power_input = compute_spectrum_power(Omega, mean_y, om_input, dom);
%     power_input_harmonics = 0;
%     for m=2:MAX_HARMONIC_N
%         power_input_harmonics = power_input_harmonics + compute_spectrum_power(Omega, mean_y, om_input * m, dom);
%     end
%     if power_input >= MIN_HARMONICS_POWER_THRESHOLD * power_input_harmonics
%         power_input = power_input + power_input_harmonics;
%     end
% 
%     Q_mean = power_input / power_natural;
%     W_mean = power_input / power_total;
% end

Omega = mean_omega;
W = zeros(Ntrials, 1);
for n=1:Ntrials
    W(n) = compute_entrainment_score(Omega, y(n, :), input_period, S);
end
W_mean = compute_entrainment_score(Omega, mean_y, input_period, S);

% mean(Q)
disp(['maximum individual entrainment score: ', num2str(max(W))]);
disp(['average individual entrainment score: ', num2str(mean(W)), ' +- ', num2str(std(W))]);
avg_upper_50_score = sort(W);
avg_upper_50_score = mean(avg_upper_50_score(round(length(avg_upper_50_score) / 2):end));
disp(['average individual upper 50% percentile score: ', num2str(avg_upper_50_score)]);
% Q_mean
disp(['complex average entrainment score: ', num2str(W_mean)]);
