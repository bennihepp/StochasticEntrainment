natural_period = 1/0.1065;
PERIOD_DEVIATION_THRESHOLD = 0.01 * natural_period;
PERIODICITY_THRESHOLD = 0.05;
PERIOD_MULTIPLE_THRESHOLD = 0.01;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
MIN_HARMONICS_POWER_THRESHOLD = 0.0;
% MAX_HARMONIC_N = 4;
MAX_HARMONIC_N = double(intmax());
entrainment_ratios = 1;

% volume = inf;
% volume = 1e6;
% volume = 1e4;
volume = 5e3; % good results for average entrainment
% volume = 1e3;
% volume = 1e2;

if volume == inf
    Ntrials = 1;
else
    Ntrials = 100;
end

dt = 1e-1;
recordStep = dt;

disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);


t0 = 0;
tf = 10000;
to = (tf - t0) / 5;
% to = tf - 50;

% % input_amplitude = 1.0;
% % input_amplitude = 0.8;
% input_amplitude = 0.0;
% input_period = 5;

% border of arnold tongue (period 15, amplitude 0.4674004, inclusive)

input_offset = 0.0;

input_period = 12;
input_amplitude = 0.2;
% input_amplitude = 0.4;

input_period = 15;
input_amplitude = 0.25;

input_period = natural_period;
% input_amplitude = 0.0;
input_amplitude = 0.1;
initial_phase = 0 * pi;
initial_phase = 1 * pi;

min_frequency = 0.01;
max_frequency = 1.0;
min_frequency = 0.0;
max_frequency = inf;

additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
multiplicative_forcing_func = @(t, x) 0;

%% simulate
tic;
% [T, output] = VanDerPol_Run(Ntrials, t0, tf, dt, volume, additive_forcing_func, multiplicative_forcing_func);
[T, output] = VanDerPol_Run2(Ntrials, t0, tf, dt, recordStep, volume, input_offset, input_amplitude, input_period, initial_phase);
forcing = input_offset + input_amplitude * sin(2*pi*T./input_period + initial_phase);
toc

% filename = ['output/simulation_Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ',_amplitude=', num2str(input_amplitude), ',_period=', num2str(input_period), '.mat'];
% save(filename);


%% plot trajectories
% figure();
% plot(T, output(:, 1));
% title(['y(1) single traces: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
% xlabel('time t');
% ylabel('state y(1)');

% figure();
% plot(T, mean(output, 2), 'LineWidth', 2.0);
% title(['y(1) average trace: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
% xlabel('time t');
% ylabel('state y(1)');

return;

%% cutoff transients
offset = find(T >= 500*input_period, 1) - 1;
% offset = 1;
T = T(offset:end-1);
output = output(offset:end-1, :);
forcing = forcing(offset:end-1);



%% substract mean
output = output - repmat(mean(output, 1), [size(output, 1), 1]);
forcing = forcing - mean(forcing);

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
[~, ff_spectrum] = compute_normalized_fft_truncated(forcing, recordStep, 2*pi*min_frequency, 2*pi*max_frequency);


%% plot phase distribution of natural mode and input mode
NUM_OF_BINS = 40;
bins = linspace(0, 2, NUM_OF_BINS);

width = 10;
height = 6;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);

% subplot(2, 1, 1);
% hold on;
% [~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ natural_period));
% hist(angle(ff_spectrum(:, ind)), bins);
% xlim([-pi, pi]);
% p = findobj(gca, 'Type', 'patch');
% set(p, 'FaceColor', 'blue', 'EdgeColor', 'black');
% hold off;
% %xlabel('phase');
% set(gca(), 'xtick', []);
% ylabel('occurence');
% %save_plot('../paper/figures/vanderpol_phase_dist_natural', h, width, height);

subplot(2, 1, 1);
xaxis_offset = 0.05;
pos = get(gca, 'Position');
set(gca, 'Position', pos + [0, +xaxis_offset, 0, -xaxis_offset]);
hold on;
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
[nelements, centers] = hist(mod(angle(y5(:, ind)), 2*pi) / (pi), bins);
bar(centers, nelements / sum(nelements), 'hist');
xlim([0, 2]);
p = findobj(gca, 'Type', 'patch');
set(p, 'FaceColor', 'blue', 'EdgeColor', 'black');
hold off;
ylabel('occurence');
set(gca(), 'xtick', []);

hx = xlabel('phase');
ylabel('occurence');

format_ticks(gca, ...
    {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'}, ...
    [], ...
    [0, 0.5, 1, 1.5, 2] ...
);
% pos = get(gca, 'Position');
% set(gca, 'Position', pos + [0, 0.05, 0, 0]);
set(hx, 'Units', 'Normalized');
pos = get(hx, 'Position');
set(hx, 'Position', pos + [0, -0.15, 0]);

pos = get(gca, 'Position');
inset = axes('Position', [pos(1) + 0.43, pos(2) + 0.11, 0.3, 0.15]);
mean_phase = mean(angle(y5(:, ind)));
window = 0.1;
inset_bins = linspace(mean_phase-window/2, mean_phase+window/2, 30);
[nelements, centers] = hist(angle(y5(:, ind)), inset_bins);
bar(centers, nelements / sum(nelements), 'hist');
xlim([mean_phase-window/2,mean_phase+window/2]);

subplot(2, 1, 2);
pos = get(gca, 'Position');
set(gca, 'Position', pos + [0, 2*xaxis_offset, 0, -xaxis_offset]);
hold on;
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
[nelements, centers] = hist(mod(angle(y(:, ind)),2*pi) / (pi), bins);
bar(centers, nelements / sum(nelements), 'hist');
xlim([0, 2]);
p = findobj(gca, 'Type', 'patch');
set(p, 'FaceColor', 'blue', 'EdgeColor', 'black');
hold off;

hx = xlabel('phase');
ylabel('occurence');

format_ticks(gca, ...
    {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'}, ...
    [], ...
    [0, 0.5, 1, 1.5, 2] ...
);
% pos = get(gca, 'Position');
% set(gca, 'Position', pos + [0, 0.00, 0, 0]);
set(hx, 'Units', 'Normalized');
pos = get(hx, 'Position');
set(hx, 'Position', pos + [0, -0.15, 0]);

pos = get(gca, 'Position');
inset = axes('Position', [pos(1) + 0.33, pos(2) + 0.11, 0.3, 0.15]);
mean_phase = mean(mod(angle(y(:, ind)), 2*pi));
window = 0.1;
inset_bins = linspace(mean_phase-window/2, mean_phase+window/2, 30);
[nelements, centers] = hist(mod(angle(y(:, ind)), 2*pi), inset_bins);
bar(centers, nelements / sum(nelements), 'hist');
xlim([mean_phase-window/2,mean_phase+window/2]);

%save_plot('../paper/figures/vanderpol_phase_dist_input', h, width, height);
save_plot([export_eps_prefix(), 'vanderpol_forcing_phase_synchronization_dist'], h, width, height);


[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
figure();
set(gca(), 'FontSize', 20);
hist(angle(ff_spectrum(:, ind)), bins);
 title(['phase distribution of forcing']);
xlabel('phase');
ylabel('occurence');

[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ natural_period));
figure();
set(gca(), 'FontSize', 20);
hist(angle(y(:, ind)), bins);
% title(['phase distribution of natural mode for volume=', num2str(volume), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');

[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
figure();
set(gca(), 'FontSize', 20);
hist(angle(y(:, ind)), bins);
% title(['phase distribution of input mode for volume=', num2str(volume), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');

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
save_plot([export_eps_prefix(), 'vanderpol_phase_dist'], h, width, height);

width = 10;
height = 4;
fontSize = 20;
h = prepare_plot(width, height, fontSize);
plot(mean_omega ./ (2 * pi), mean(abs(y), 1) .^ 2, 'Color', 'red', 'LineWidth', 1.5);
xlabel('frequency f');
ylabel('power |y|^2');
save_plot('../paper/figures/vanderpol_single_spectrum', h, width, height);


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
% 
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
% 
% end

S = struct();
S.natural_period = natural_period;
S.FREQUENCY_NEIGHBOURHOOD_FACTOR = FREQUENCY_NEIGHBOURHOOD_FACTOR;
S.MAX_HARMONIC_N = MAX_HARMONIC_N;
S.MIN_HARMONICS_POWER_THRESHOLD = MIN_HARMONICS_POWER_THRESHOLD;
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
