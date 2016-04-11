%% some functions in the parent folder are used
addpath('../');
addpath('../plotting/');


% not necessary for deterministic simulation
% omega = 600;
% Ntrials = 100;

t0 = 0.0;
tf = 400*24;
recordStep = (tf - t0) / 5000.0;
to = (tf - t0) / 10;


input_offset = 1.0;
input_period = 30.0;
input_amplitude = 0.28;

initial_phases = linspace(0, 2*pi, 20);

min_frequency = 0.0;
max_frequency = inf;


%% simulate

output_phases = zeros(size(initial_phases));
display('loop');
for i=1:length(initial_phases)
    display(['  i=', int2str(i)]);
    initial_phase = initial_phases(i);

    [T, output] = RunODE(t0, tf, recordStep, ...
        input_offset, input_amplitude, input_period, initial_phase);
%     [T, output, ~] = Run(Ntrials, t0, tf, recordStep, omega, ...
%         input_offset, input_amplitude, input_period, initial_phase);

    output = output - repmat(mean(output, 1), [size(output, 1), 1]);

    %% compute spectras
    [omega, ff_spectrum] = compute_normalized_fft_truncated(output, recordStep, 2*pi*min_frequency, 2*pi*max_frequency);

    [~, ind] = min(abs(omega ./ (2 * pi) - 1 ./ input_period));
    output_phase = angle(ff_spectrum(ind));
    output_phases(i) = output_phase;

end

%% plot mapping of input and output phase

X = initial_phases;
Y = output_phases;
Y(Y < Y(1)) = Y(Y < Y(1)) + 2*pi;
% Y(end) = Y(end) + 2*pi;
X = X / (pi);
Y = (Y - Y(1)) / (pi);

xaxis_offset = 0.05;
yaxis_offset = 0.08;
width = 6;
height = 4;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
hold on;
plot(X, Y, 'b', 'LineWidth', 2.0);
hold off;
hx = xlabel('input phase \phi');
hy = ylabel('output phase \theta');
set(gca, 'xlim', [0, 2], 'ylim', [0, 2]);
pos = get(gca, 'Position');
set(gca, 'Position', pos + [yaxis_offset, 0, -yaxis_offset, 0]);

format_ticks(gca, ...
    {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'}, ...
    {'\theta_{0}+0', '\theta_{0}+\pi', '\theta_{0}+2\pi'}, ...
    [0, 0.5, 1, 1.5, 2], ...
    [0, 1, 2]...
);

set(hx, 'Units', 'Normalized');
pos = get(hx, 'Position');
set(hx, 'Position', pos + [0, -xaxis_offset, 0]);

set(hy, 'Units', 'Normalized');
pos = get(hy, 'Position');
set(hy, 'Position', pos + [-yaxis_offset-0.05, 0, 0]);

% save_plot([export_eps_prefix(), 'leloup_goldbeter_circadian_forcing_phase_synchronization'], h, width, height);
