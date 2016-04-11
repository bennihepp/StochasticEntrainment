%% parameters for scoring entrainment
natural_period = 23.8607;
PERIOD_DEVIATION_THRESHOLD = 0.01 * natural_period;
PERIODICITY_THRESHOLD = 0.05;
PERIOD_MULTIPLE_THRESHOLD = 0.01;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;
MIN_HARMONICS_POWER_THRESHOLD = 0.0;
% MAX_HARMONIC_N = 4;
MAX_HARMONIC_N = double(intmax());
entrainment_ratios = 1;


%% parameters for the simulation

% scaling constant for the system
omega = 600;

% number of trajectories to simulate (for infinite volume only one trajectory is necessary)
Ntrials = 100;


% initial time
t0 = 0.0;
% final time
tf = 400*24;
% offset time to cutoff to reduce transient effects
to = (tf - t0) / 10;


% time-interval for saving of the output state
recordStep = (tf - t0)/5000;

disp([' Ntrials=', int2str(Ntrials), ' recordStep=', num2str(recordStep)]);


%% parameters for the forcing function

input_offset = 1.0;
% input_period = 24.0;
% input_amplitude = 0.2;
input_period = 30.0;
input_amplitude = 0.28;


%% list of initial phases to compare (the plotting code only works with 2 phases)
initial_phases = [0 * pi, 1 * pi];
number_of_initial_phases = length(initial_phases);


%% parameters for computation of spectra

% minimum and maximum frequency to consider in the fourier spectrum
min_frequency = 0.0;
max_frequency = inf;


%% simulate

all_Ts = cell(number_of_initial_phases, 1);
all_outputs = cell(number_of_initial_phases, 1);

for i=1:number_of_initial_phases
    initial_phase = initial_phases(i);

    % run MATLAB model and solver
%     [T, output] = VanDerPol_Run(Ntrials, t0, tf, dt, volume, additive_forcing_func, multiplicative_forcing_func);

    % run Java model and solver
    printMessages = true;
    [T, output, ~] = Run(Ntrials, t0, tf, recordStep, omega, ...
        input_offset, input_amplitude, input_period, initial_phase, printMessages);
    all_Ts{i} = T;
    all_outputs{i} = output;

%     %% plot trajectories
%     figure();
%     plot(T, output(:, 1));
%     title(['y(1) single traces: Ntrials=', int2str(Ntrials), ' recordStep=', num2str(recordStep), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
%     xlabel('time t');
%     ylabel('state y(1)');
% 
%     figure();
%     plot(T, mean(output, 2), 'LineWidth', 2.0);
%     title(['y(1) average trace: Ntrials=', int2str(Ntrials), ' recordStep=', num2str(recordStep), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
%     xlabel('time t');
%     ylabel('state y(1)');

end


%% save unprocessed data
filename = ['output/phase_simulations_Ntrials=', int2str(Ntrials), ...
    '_recordStep=', num2str(recordStep), '_omega=', num2str(omega), ...
    '_offset=', num2str(input_offset), ...
    '_amplitude=', num2str(input_amplitude), ...
    '_period=', num2str(input_period), '.mat'];
save(filename);

return;


%% some functions in the parent folder are used
addpath('../');
addpath('../plotting/');


%% cutoff transients
for i=1:number_of_initial_phases
    T = all_Ts{i};
    output = all_outputs{i};

    offset = find(T >= 20*input_period, 1) - 1;
    T = T(offset:end-1);
    output = output(offset:end-1, :);

    all_Ts{i} = T;
    all_outputs{i} = output;
end


%% substract mean

for i=1:number_of_initial_phases
    output = all_outputs{i};

    output = output - repmat(mean(output, 1), [size(output, 1), 1]);

    all_outputs{i} = output;
end


%% compute spectra

all_omegas = cell(number_of_initial_phases, 1);
all_ys = cell(number_of_initial_phases, 1);
all_mean_omegas = cell(number_of_initial_phases, 1);

for i=1:number_of_initial_phases
    output = all_outputs{i};

    omega = [];
    y = [];
    for j=Ntrials:-1:1
        [omega1, y1] = compute_normalized_fft_truncated(output(:,j)', recordStep, 2*pi*min_frequency, 2*pi*max_frequency);
        omega = [omega; omega1];
        y = [y; y1];
    end

%     mean_y = mean(y, 1);
    mean_omega = mean(omega, 1);

    all_omegas{i} = omega;
    all_ys{i} = y;
    all_mean_omegas{i} = mean_omega;
end


%% plot phase distribution of the forcing frequency mode
%% with the first and second initial input phase

NUM_OF_BINS = 40;
bins = linspace(0, 2, NUM_OF_BINS);

width = 10;
height = 6;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);

mean_omega = all_mean_omegas{1};
y = all_ys{1};
subplot(2, 1, 1);
xaxis_offset = 0.05;
pos = get(gca, 'Position');
set(gca, 'Position', pos + [0, +xaxis_offset, 0, -xaxis_offset]);
hold on;
[~, ind] = min(abs(mean_omega ./ (2 * pi) - 1 ./ input_period));
[nelements, centers] = hist(mod(angle(y(:, ind)), 2*pi) / (pi), bins);
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
inset = axes('Position', [pos(1) + 0.3, pos(2) + 0.11, 0.3, 0.15]);
mean_phase = mean(angle(y(:, ind)));
window = 0.1;
inset_bins = linspace(mean_phase-window/2, mean_phase+window/2, 30);
[nelements, centers] = hist(angle(y(:, ind)), inset_bins);
bar(centers, nelements / sum(nelements), 'hist');
xlim([mean_phase-window/2,mean_phase+window/2]);


mean_omega = all_mean_omegas{2};
y = all_ys{2};
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
inset = axes('Position', [pos(1) + 0.45, pos(2) + 0.11, 0.3, 0.15]);
mean_phase = mean(mod(angle(y(:, ind)), 2*pi));
mean_phase = mean(angle(y(:, ind)));
window = 0.1;
inset_bins = linspace(mean_phase-window/2, mean_phase+window/2, 30);
[nelements, centers] = hist(mod(angle(y(:, ind)), 2*pi), inset_bins);
bar(centers, nelements / sum(nelements), 'hist');
xlim([mean_phase-window/2,mean_phase+window/2]);

% save_plot([export_eps_prefix(), 'leloup_goldbeter_circadian_forcing_phase_synchronization_dist'], h, width, height);
