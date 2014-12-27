natural_period = 23.7473;
PERIOD_DEVIATION_THRESHOLD = 0.01 * natural_period;
PERIODICITY_THRESHOLD = 0.05;
PERIOD_MULTIPLE_THRESHOLD = 0.01;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
MIN_HARMONICS_POWER_THRESHOLD = 0.0;
MAX_HARMONIC_N = 4;
entrainment_ratios = 1:2;
addpath([getenv('HOME'), '/Documents/MATLAB/plotting']);

volume = inf;
% volume = 1e-20;
% volume = 2e-18;

if volume == inf
    Ntrials = 1;
else
    Ntrials = 1000;
    Ntrials = 3;
end

dt = 0.002;
recordStep = 400 * dt;

disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);

t0 = 0;
tf = 200*72;
to = (tf - t0) / 5;

input_offset = 1.0;

% input_period = 36.0;
% input_amplitude = 0.5;

input_period = 30.0;
input_amplitude = 0.18;

% input_period = 24.0;
% input_amplitude = 0.01;

% input_period = 35.5;
% input_amplitude = 0.0;


%% simulate
tic;
printMessages = true;
[T, output] = CircadianClock_Run(Ntrials, t0, tf, dt, recordStep, volume, input_offset, input_amplitude, input_period, printMessages);
toc

% filename = ['output/simulation_Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), '_offset=', num2str(input_offset), ',_amplitude=', num2str(input_amplitude), ',_period=', num2str(input_period), '.mat'];
% save(filename);

return;

%% plot trajectories
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

%% cutoff transients
% offset_time = (tf - t0) / 5;
% offset_time = max(offset_time, 1000);
offset_time = to;
offset = find(T >= offset_time, 1);
T = T(offset:end);
output = output(offset:end, :);


%% plot traces after transients

w = find(T > T(end) - 400, 1);
q = find(T > T(w) - 200, 1);
TT = T(q:w);
TT = TT - TT(1);
trunc = output(q:w, :);

% figure;
% %     plot(T, output(:, i), 'Color', cmap(i, :), 'LineWidth', 2.0);
% plot(TT, trunc(:, 1), 'LineWidth', 2.0);
% title(['y(1) single trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
% xlabel('time t');
% ylabel('state y');

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

width = 10;
height = 3;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
hold on;
%cmap = colormap('Lines');
color1 = [0, 0, 1.0];
color2 = [0, 1.0, 1.0];
color3 = [0, 1.0, 0.0];
% average_color = [1.0, 0.5, 0.0];
average_color = [241, 140, 22] / 255;
cmap = [color1; color2; color3];
for i=1:3
    plot(TT, trunc(:, i), 'Color', cmap(i, :), 'LineWidth', 0.5);
end
plot(TT, mean(trunc, 2), '-', 'Color', average_color, 'LineWidth', 2.0);
xlabel('time t');
ylabel('state y');
hold off;
save_plot([export_eps_prefix(), 'circadian_average_and_single_trace'], h, width, height);

% width = 10;
% height = 3;
% fontSize = 0.5 * (width * height);
% h = prepare_plot(width, height, fontSize);
% hold on;
% deterministic_color = [0.9, 0.0, 0.0];
% plot(TT, mean(trunc, 2), '-', 'Color', deterministic_color, 'LineWidth', 1.0);
% xlabel('time t');
% ylabel('state y');
% hold off;
% save_plot([export_eps_prefix(), 'vanderpol_deterministic_trace'], h, width, height);


%% substract mean
output = output - repmat(mean(output, 1), [size(output, 1), 1]);

%% compute spectras
addpath('../');

% min_frequency = 0.005;
% max_frequency = 0.5;
min_frequency = 0.0;
max_frequency = inf;

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
title(['y(1) first trace fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('frequency f');
ylabel('power |y|^2');

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

NUM_OF_BINS = 100;
bins = linspace(-pi, pi, NUM_OF_BINS);
width = 10;
height = 4;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
subplot(2, 1, 1);
hold on;
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ natural_period));
hist(angle(y(:, ind)), bins);
p = findobj(gca, 'Type', 'patch');
set(p, 'FaceColor', 'blue', 'EdgeColor', 'black');
hold off;
%xlabel('phase');
set(gca(), 'xtick', []);
ylabel('occurence');
%save_plot('../paper/figures/vanderpol_phase_dist_natural', h, width, height);
subplot(2, 1, 2);
hold on;
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
hist(angle(y(:, ind)), bins);
p = findobj(gca, 'Type', 'patch');
set(p, 'FaceColor', 'blue', 'EdgeColor', 'black');
hold off;
xlabel('phase');
ylabel('occurence');
%save_plot('../paper/figures/vanderpol_phase_dist_input', h, width, height);
save_plot([export_eps_prefix(), 'circadian_phase_dist'], h, width, height);


%% compute autocorrelation if only one trajectory is simulated
if Ntrials == 1
    corr = xcorr(output - mean(output, 2), 'unbiased');
    figure();
    plot(corr);
    title('autocorrelation');
    [pks, locs] = findpeaks(corr);
    peak_distances = locs(2:end) - locs(1:end-1);
    mean_peak_distance = mean(peak_distances);
    std_peak_distance = std(peak_distances);

    mean_period = mean_peak_distance * dt;
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

S = struct();
S.natural_period = natural_period;
S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
S.MAX_HARMONIC_N = MAX_HARMONIC_N;
S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
S.MIN_HARMONICS_POWER_THRESHOLD = 0;
S.entrainment_ratios = entrainment_ratios;

Omega = mean_omega;
W = zeros(Ntrials, 1);
for n=1:Ntrials
    W(n) = compute_entrainment_score(Omega, y(n, :), input_period, S);
end
W_mean = compute_entrainment_score(Omega, mean_y, input_period, S);

% mean(Q)
disp(['maximum individual entrainment score: ', num2str(max(W))]);
disp(['average individual entrainment score: ', num2str(mean(W)), ' +- ', num2str(std(W))]);
% Q_mean
disp(['complex average entrainment score: ', num2str(W_mean)]);
