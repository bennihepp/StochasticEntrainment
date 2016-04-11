%% some functions in the parent folder are used
addpath('../');
addpath('../plotting/');


%% parameters for scoring entrainment
natural_period = 23.8607;
PERIOD_DEVIATION_THRESHOLD = 0.01 * natural_period;
PERIODICITY_THRESHOLD = 0.05;
PERIOD_MULTIPLE_THRESHOLD = 0.01;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
MIN_HARMONICS_POWER_THRESHOLD = 0.0;
MAX_HARMONIC_N = 4;
entrainment_ratios = 1:2;


%% parameters for the simulation

% scaling constant for the system
omega = 600;

% number of trajectories to simulate (for infinite volume only one trajectory is necessary)
Ntrials = 4;
% Ntrials = 8;
Ntrials = 500;
% Ntrials = 3;

% time-interval for saving of the output state
% recordStep = (tf - t0)/5000;
recordStep = (tf - t0)/10000;

disp([' Ntrials=', int2str(Ntrials), ' recordStep=', num2str(recordStep)]);


% initial time
t0 = 0;
% final time
%tf = 100 * 24;
tf = 50 * 24;
% tf = 400 * 24;
% offset time to cutoff to reduce transient effects
to = (tf - t0) / 10;


%% parameters for the forcing function

input_offset = 1.0;
initial_phase = 0.0;

% % entrained
% input_period = 30.0;
% input_amplitude = 0.3;

% % not entrained
% input_period = 30.0;
% input_amplitude = 0.2;

input_period = 28;
input_amplitude = 0.2;


%% parameters for computation of spectra

% minimum and maximum frequency to consider in the fourier spectrum
min_frequency = 0.0;
max_frequency = 1 / 3;
% min_frequency = 0.0;
% max_frequency = inf;


%% initialize options structure
S = struct();
S.natural_period = natural_period;
S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
S.MAX_HARMONIC_N = MAX_HARMONIC_N;
S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
S.MIN_HARMONICS_POWER_THRESHOLD = 0;
S.entrainment_ratios = entrainment_ratios;


%% simulate
tic;
printMessages = true;
[T, output, ~] = Run(Ntrials, t0, tf, recordStep, omega, ...
    input_offset, input_amplitude, input_period, initial_phase, printMessages);
toc

% filename = ['output/simulation_Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), '_offset=', num2str(input_offset), ',_amplitude=', num2str(input_amplitude), ',_period=', num2str(input_period), '.mat'];
% save(filename);

return;

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


%% plot traces after transients

w = find(T > T(end) - 400, 1);
q = find(T > T(w) - 150, 1);
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
title(['y(1) single traces: Ntrials=', int2str(Ntrials)]);
xlabel('time t');
ylabel('state y');

figure();
% plot(T, mean(output, 2), 'LineWidth', 2.0);
plot(TT, mean(trunc, 2), 'LineWidth', 2.0);
title(['y(1) average trace: Ntrials=', int2str(Ntrials), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('time t');
ylabel('state y(1)');


width = 10;
height = 5;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
%cmap = colormap('Lines');
color1 = [0, 0, 1.0];
color2 = [0, 1.0, 1.0];
color3 = [0, 1.0, 0.0];
% average_color = [1.0, 0.5, 0.0];
average_color = [241, 140, 22] / 255;
cmap = [color1; color2; color3];
subplot(2, 1, 1);
hold on;
for i=1:3
    plot(TT, trunc(:, i) / omega, 'Color', cmap(i, :), 'LineWidth', 0.5);
end
plot(TT, mean(trunc, 2) / omega, '-', 'Color', average_color, 'LineWidth', 2.0);
hold off;
ylabel('state y');
set(gca(), 'xtick', []);
box off;
subplot(2, 1, 2);
deterministic_color = [0.9, 0.0, 0.0];
plot(TT_det, mean(trunc_det, 2), '-', 'Color', deterministic_color, 'LineWidth', 1.0);
xlabel('time t');
ylabel('state y');
box off;
% save_plot([export_eps_prefix(), 'leloup_goldbeter_circadian_average_and_single_trace'], h, width, height);

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

% compute the spectrum for each trajectory
omega = [];
y = [];
for i=1:Ntrials
    display([int2str(i), ' out of ', int2str(Ntrials)]);
    [omega1, y1] = compute_normalized_fft_truncated(output(:,i)', recordStep, 2*pi*min_frequency, 2*pi*max_frequency);
    omega = [omega; omega1];
    y = [y; y1];
end

% compute the average spectrum (i.e. the spectrum of the average)
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
% save_plot([export_eps_prefix(), 'leloup_goldbeter_circadian_phase_dist'], h, width, height);


%% compute entrainment scores

% compute score for each trajectory
W = zeros(Ntrials, 1);
for n=1:Ntrials
    W(n) = compute_entrainment_score(mean_omega, y(n, :), input_period, S);
end
% compute score of average trajectory
W_mean = compute_entrainment_score(mean_omega, mean_y, input_period, S);

disp(['maximum individual entrainment score: ', num2str(max(W))]);
disp(['average individual entrainment score: ', num2str(mean(W)), ' +- ', num2str(std(W))]);

% compute the average score of the upper 50 percentile
W_sort = sort(W);
W_upper_50 = W_sort(round(length(W_sort) / 2):end);
W_avg_upper_50 = mean(W_upper_50);
disp(['individual entrainment score: ', num2str(W_avg_upper_50)]);

disp(['complex average entrainment score: ', num2str(W_mean)]);
