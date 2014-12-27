addpath([getenv('HOME'), '/Documents/MATLAB/plotting/']);


A_0 = 1.0;
A_1 = 1.0;

omega_0 = 3;
omega_1 = 7;
epsilon = 0.3;

t0 = 0;
% tf = 3;
tf = 500;
to = (tf - t0) / 5;
% numOfSteps = 1000;
% dt = tf / numOfSteps;
dt = 0.01;
Ntrials = 1000;

disp(['Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);

% phi_0_mean = 1.0 * randn(1, Ntrials) * pi;
phi_0_mean = 1.0 * (2 * rand(1, Ntrials) - 1.0) * pi;
phi_0_std = 0.0 * pi;
phi_0_rand = @(n, s1, s2) phi_0_mean(n) + phi_0_std * randn(s1, s2);
phi_0_det = @(n, s1, s2) 0;
% phi_0_rand = @(n1, n2) (2 * rand(n1, n2) - 1) * pi;
% phi_1_mean = 0.1 * randn(1, Ntrials) * pi;
phi_1_mean = 0.5 * (2 * rand(1, Ntrials) - 1.0) * pi;
phi_1_std = 0.0 * pi;
phi_1_rand = @(n, s1, s2) phi_1_mean(n) + phi_1_std * randn(s1, s2);
phi_1_det = @(n, s1, s2) 0;


%% simulate
tic;
[T_det, output_det] = Heuristic_Run(1, t0, tf, dt, A_0, A_1, epsilon, omega_0, omega_1, phi_0_det, phi_1_det);
[T, output] = Heuristic_Run(Ntrials, t0, tf, dt, A_0, A_1, epsilon, omega_0, omega_1, phi_0_rand, phi_1_rand);
toc


%% cutoff transients
% offset_time = (tf - t0) / 5;
% offset_time = max(offset_time, 1000);
offset_time = to;
offset = find(T >= offset_time, 1);
T = T(offset:end);
output = output(offset:end, :);
T_det = T_det(offset:end);
output_det = output_det(offset:end, :);


%% only plot a short amount of time
TT = T - T(1);
TT = T(1:1000);
output_tmp = output(1:1000, :);
TT_det = T_det - T_det(1);
TT_det = T_det(1:1000);
output_det_tmp = output_det(1:1000, :);


%% plot trajectories
figure();
hold on;
cmap = colormap('Lines');
for i=1:3
    plot(TT, output_tmp(:, i), 'Color', cmap(i, :), 'LineWidth', 1.0);
end
%title(['y single traces: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
xlabel('time t');
ylabel('state y');
plot(TT, mean(output_tmp, 2), '-', 'Color', 0.3*ones(3,1), 'LineWidth', 4.0);
hold off;

% width = 10;
% height = 3;
% fontSize = 0.5 * (width * height);
% h = prepare_plot(width, height, fontSize);
% hold on;
% % cmap = colormap('Lines');
% color1 = [0, 0, 1.0];
% color2 = [0, 1.0, 1.0];
% color3 = [0, 1.0, 0.0];
% % average_color = [1.0, 0.5, 0.0];
% average_color = [241, 140, 22] / 255;
% cmap = [color1; color2; color3];
% for i=1:3
%     plot(TT, output_tmp(:, i), 'Color', cmap(i, :), 'LineWidth', 1.0);
% end
% plot(TT, mean(output_tmp, 2), '-', 'Color', average_color, 'LineWidth', 2.0);
% xlabel('time t');
% ylabel('state y');
% hold off;
% save_plot([export_eps_prefix(), 'heuristic_average_and_single_trace'], h, width, height);

width = 10;
height = 3;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
hold on;
% cmap = colormap('Lines');
color4 = [1.0, 0, 1.0];
plot(TT, output_det_tmp(:, 1), 'Color', color4, 'LineWidth', 1.0);
xlabel('time t');
ylabel('state y');
hold off;
save_plot([export_eps_prefix(), 'heuristic_deterministic_trace'], h, width, height);

width = 10;
height = 3;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
hold on;
% cmap = colormap('Lines');
color1 = [0, 0, 1.0];
color2 = [0, 1.0, 1.0];
color3 = [0, 1.0, 0.0];
% average_color = [1.0, 0.5, 0.0];
average_color = [241, 140, 22] / 255;
cmap = [color1; color2; color3];
for i=1:3
    plot(TT, output_tmp(:, i), 'Color', cmap(i, :), 'LineWidth', 1.0);
end
xlabel('time t');
ylabel('state y');
hold off;
save_plot([export_eps_prefix(), 'heuristic_single_trace'], h, width, height);

width = 10;
height = 3;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
hold on;
average_color = [241, 140, 22] / 255;
plot(TT, mean(output_tmp, 2), '-', 'Color', average_color, 'LineWidth', 1.0);
xlabel('time t');
ylabel('state y');
hold off;
ylim([-2, 2]);
save_plot([export_eps_prefix(), 'heuristic_average_trace'], h, width, height);


%% compute spectras
addpath('../');

min_frequency = 0.0;
% max_frequency = 0.19;
max_frequency = 1.99;
% max_frequency = inf;

omega = [];
y = [];
for i=Ntrials:-1:1
    i
    [omega1, y1] = compute_normalized_fft_truncated(output(:,i)', dt, 2*pi*min_frequency, 2*pi*max_frequency);
    omega = [omega; omega1];
    y = [y; y1];
end
mean_y = mean(y, 1);
mean_omega = mean(omega, 1);
[omega_dt, y_det] = compute_normalized_fft_truncated(output_det(:,1)', dt, 2*pi*min_frequency, 2*pi*max_frequency);

%% plot spectras
figure();
hold on;
cmap = colormap('Lines');
for i=1:3
    plot(mean_omega ./ (2 * pi), mean(abs(y(i,:)), 1) .^ 2, 'Color', cmap(i, :), 'LineWidth', 2.0);
end
hold off;
title(['y single traces fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
xlabel('time t');
ylabel('state y');

figure();
plot(mean_omega ./ (2 * pi), mean(abs(y(1,:)), 1) .^ 2);
title(['y(1) first trace fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
xlabel('frequency f');
ylabel('power |y|^2');

figure();
plot(mean_omega ./ (2 * pi), abs(mean_y) .^ 2, 'LineWidth', 2.0);
%title(['y(1) complex average fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
xlabel('frequency f');
ylabel('power |y|^2');

figure();
plot(mean_omega ./ (2 * pi), mean(abs(y), 1) .^ 2);
%title(['y(1) absolute average fft: Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);
xlabel('frequency f');
ylabel('power |y|^2');

width = 10;
height = 6;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
values = abs(y_det) .^ 2;
values = values / sum(values);
plot(mean_omega ./ (2 * pi), values, 'Color', color4, 'LineWidth', 1.5);
xlabel('frequency f');
ylabel('power |y|^2');
% ylim([0, 0.5]);
save_plot([export_eps_prefix(), 'heuristic_deterministic_spectrum'], h, width, height);

width = 10;
height = 6;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
values = mean(abs(y), 1) .^ 2;
values = values / sum(values);
plot(mean_omega ./ (2 * pi), values, 'Color', color2, 'LineWidth', 1.5);
xlabel('frequency f');
ylabel('power |y|^2');
% ylim([0, 0.5]);
save_plot([export_eps_prefix(), 'heuristic_single_spectrum'], h, width, height);

width = 10;
height = 6;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
values = abs(mean_y) .^ 2;
values = values / sum(values);
plot(mean_omega ./ (2 * pi), values, 'Color', average_color, 'LineWidth', 1.5);
xlabel('frequency f');
ylabel('power |y|^2');
% ylim([0, 0.5]);
save_plot([export_eps_prefix(), 'heuristic_average_spectrum'], h, width, height);

return;


%% plot phase distribution of natural mode and input mode
NUM_OF_BINS = 50;
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ natural_period));
figure();
hist(angle(y(:, ind)), NUM_OF_BINS);
title(['phase distribution of natural mode for volume=', num2str(volume), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');

[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
figure();
hist(angle(y(:, ind)), NUM_OF_BINS);
title(['phase distribution of input mode for volume=', num2str(volume), ', input period=', num2str(input_period), ', input amplitude=', num2str(input_amplitude), ', Ntrials=', int2str(Ntrials)]);
xlabel('phase');
ylabel('occurence');


% filename = ['output/simulation_Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' offset=', num2str(TNF_offset), ' amplitude=', num2str(TNF_amplitude), ' period=', num2str(TNF_period), '.mat'];
% save(filename);

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
Omega = mean_omega;
Q = zeros(Ntrials, 1);
W = zeros(Ntrials, 1);
om_natural = 2 * pi / natural_period;
om_input = 2 * pi / input_period;
dom = FREQUENCY_NEIGHBOURHOOD_FACTOR * om_natural;
if abs(om_natural - om_input) < 2 * dom
    Q(:) = inf;
    W(:) = inf;
else
    for n=1:Ntrials
        power_total = sum(abs(y(n,:)).^2);
        power_natural = compute_spectrum_power(Omega, y(n,:), om_natural, dom);
%         for n=2:MAX_HARMONIC_N
%             power_natural = power_natural + compute_spectrum_power(Omega, y(n,:), om_natural * n, dom);
%         end
        power_input = compute_spectrum_power(Omega, y(n,:), om_input, dom);
        power_input_harmonics = 0;
        for m=2:MAX_HARMONIC_N
            power_input_harmonics = power_input_harmonics + compute_spectrum_power(Omega, y(n,:), om_input * m, dom);
        end
        if power_input >= MIN_HARMONICS_POWER_THRESHOLD * power_input_harmonics
            power_input = power_input + power_input_harmonics;
        end

        Q(n) = power_input / power_natural;
        W(n) = power_input / power_total;

    end

end

if abs(om_natural - om_input) < 2 * dom
    Q_mean = inf;
    W_mean = inf;
else
    power_total = sum(abs(mean_y).^2);
    power_natural = compute_spectrum_power(Omega, mean_y, om_natural, dom);
%     for n=2:MAX_HARMONIC_N
%         power_natural = power_natural + compute_spectrum_power(Omega, mean_y, om_natural * n, dom);
%     end
    power_input = compute_spectrum_power(Omega, mean_y, om_input, dom);
    power_input_harmonics = 0;
    for m=2:MAX_HARMONIC_N
        power_input_harmonics = power_input_harmonics + compute_spectrum_power(Omega, mean_y, om_input * m, dom);
    end
    if power_input >= MIN_HARMONICS_POWER_THRESHOLD * power_input_harmonics
        power_input = power_input + power_input_harmonics;
    end

    Q_mean = power_input / power_natural;
    W_mean = power_input / power_total;

end

% mean(Q)
disp(['maximum individual entrainment score: ', num2str(max(W))]);
disp(['average individual entrainment score: ', num2str(mean(W)), ' +- ', num2str(std(W))]);
% Q_mean
disp(['complex average entrainment score: ', num2str(W_mean)]);
