natural_period = 1/0.1065;
PERIOD_DEVIATION_THRESHOLD = 0.01 * natural_period;
PERIODICITY_THRESHOLD = 0.05;
PERIOD_MULTIPLE_THRESHOLD = 0.01;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
MIN_HARMONICS_POWER_THRESHOLD = 1.0;
MAX_HARMONIC_N = 4;

volume = inf;
% volume = 1e6;
% volume = 1e4;
% volume = 5e3; % good results for average entrainment
% volume = 1e3;
% volume = 1e2;

if volume == inf
    Ntrials = 1;
    dt = 1e-1;
else
    Ntrials = 1000;
    dt = 1e-1;
end

disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);


t0 = 0;
tf = 10000;
to = (tf - t0) / 5;

% % input_amplitude = 1.0;
% % input_amplitude = 0.8;
% input_amplitude = 0.0;
% input_period = 5;

% border of arnold tongue (period 15, amplitude 0.4674004, inclusive)

input_period = 15;
input_amplitude = 0.3;

additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
multiplicative_forcing_func = @(t, x) 0;

%% simulate
tic;
[T, output] = VanDerPol_Run(Ntrials, t0, tf, dt, volume, additive_forcing_func, multiplicative_forcing_func);
toc

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


%% compute spectras
addpath('../');

omega = [];
y = [];
for i=Ntrials:-1:1
    min_frequency = 0.01;
    max_frequency = 1.0;
    [omega1, y1] = compute_normalized_fft_truncated(output(:,i)', dt, 2*pi*min_frequency, 2*pi*max_frequency);
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
