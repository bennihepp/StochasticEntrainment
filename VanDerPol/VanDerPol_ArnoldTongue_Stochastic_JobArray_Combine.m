output_prefix = 'output/VanDerPol_ArnoldTongue_Stochastic_JobArray';

addpath('../');

% natural_period = 1/0.1065;

volume = 5e3;

if volume == inf
    Ntrials = 1;
    dt = 1e-1;
else
    Ntrials = 100;
    dt = 1e-1;
end

disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);

t0 = 0;
tf = 10000;


input_amplitudes = 0.0:0.1:1.5;
input_periods = 1:0.5:20.0;


min_frequency = 0.01;
max_frequency = 1.0;

T = t0:dt:tf;
offset_time = (tf - t0) / 5;
offset_time = min(offset_time, 1000);
offset = find(T >= offset_time, 1);
T = T(offset:end);
[Omega, ~] = compute_normalized_fft_truncated(T, dt, 2*pi*min_frequency, 2*pi*max_frequency);

num_of_iterations = length(input_periods) * length(input_amplitudes);
per_ind = zeros(num_of_iterations, 1);
amp_ind = zeros(num_of_iterations, 1);
l = 1;
for i=1:length(input_periods)
    for j=1:length(input_amplitudes)
        per_ind(l) = i;
        amp_ind(l) = j;
        l = l + 1;
    end
end

Y = zeros(length(input_periods), length(input_amplitudes), length(Omega));
Yabs = zeros(length(input_periods), length(input_amplitudes), length(Omega));
PDmean = zeros(length(input_periods), length(input_amplitudes));
PDstd = zeros(length(input_periods), length(input_amplitudes));

for l=1:num_of_iterations
    display(['l=', int2str(l), ' out of ', int2str(num_of_iterations)]);
    i = per_ind(l);
    j = amp_ind(l);
    input_period = input_periods(i);
    input_amplitude = input_amplitudes(j);
    display(['input_period=', num2str(input_period), ' input_amplitude=', num2str(input_amplitude)]);
    filename = [output_prefix, '_l=1', int2str(l), '.mat'];
    S = load(filename);
    Y(i, j, :) = S.Y;
    Yabs(i, j, :) = S.Yabs;
    PDmean(i, j) = S.PDmean;
    PDstd(i, j) = S.PDstd;
end

S = struct();
S.natural_period = natural_period;
S.volume = volume;
S.t0 = t0;
S.tf = tf;
S.dt = dt;
S.Ntrials = Ntrials;
S.input_periods = input_periods;
S.input_amplitudes = input_amplitudes;
S.min_frequency = min_frequency;
S.max_frequency = max_frequency;
S.T = T;
S.Omega = Omega;
S.Y = Y;
S.Yabs = Yabs;
% S.X = X;
S.PDmean = PDmean;
S.PDstd = PDstd;

date_string = datestr(clock());
filename = ['VanDerPol_ArnoldTongue_Stochastic_volume=', num2str(volume), '_', date_string, '.mat'];
save(filename, '-struct', 'S');
