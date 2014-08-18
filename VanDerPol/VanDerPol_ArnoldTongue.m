addpath('../');

natural_period = 1/0.1065;

volume = inf;
omega = volume;

if volume == inf
    Ntrials = 1;
    dt = 1e-1;
end

t0 = 0;
tf = 10000;


% input_amplitudes = 0.0:0.2:2.0;
% input_periods = 5:2.0:40;

input_amplitudes = 0.0:0.05:1.5;
input_periods = 5:1.0:30;

% input_amplitudes = 0.0:0.1:1.0;
% input_periods = 1:1:20;


min_frequency = 0.01;
max_frequency = 1.0;

T = t0:dt:tf;
[Omega, ~] = compute_normalized_fft_truncated(T, dt, 2*pi*min_frequency, 2*pi*max_frequency);
% X = zeros(length(input_periods), length(input_amplitudes), length(T));
Y = zeros(length(input_periods), length(input_amplitudes), length(Omega));
PDmean = zeros(length(input_periods), length(input_amplitudes));
PDstd = zeros(length(input_periods), length(input_amplitudes));


for i=1:length(input_periods)
% parfor i=1:length(input_periods)
    display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
    input_period = input_periods(i);

%     XX = zeros(length(input_amplitudes), length(T));
    YY = zeros(length(input_amplitudes), length(Omega));
    PDmean_tmp = zeros(length(input_amplitudes), 1);
    PDstd_tmp = zeros(length(input_amplitudes), 1);

    for j=1:length(input_amplitudes)
        display(['i=', int2str(i), ' of ', int2str(length(input_periods)), ', j=', int2str(j), ' of ', int2str(length(input_amplitudes))]);
        input_amplitude = input_amplitudes(j);

        additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
        multiplicative_forcing_func = @(t, x) 0;

        [TT, output] = VanDerPol_Run(Ntrials, t0, tf, dt, omega, additive_forcing_func, multiplicative_forcing_func);

%         XX(j, :) = output;

        offset_time = (tf - t0) / 5;
        offset_time = min(offset_time, 1000);
        offset = find(TT >= offset_time, 1);
        TT = TT(offset:end);
        output = output(offset:end, :);

        %% Fourier spectrum analysis
        [omega1, y1]= compute_normalized_fft_truncated(output, dt, 2*pi*min_frequency, 2*pi*max_frequency);
%         Omega = omega;
        YY(j, :) = y1;

        %% Autocorrelation analysis
        corr = xcorr(output - mean(output), 'unbiased');
        [pks, locs] = findpeaks(corr);
        peak_distances = locs(2:end) - locs(1:end-1);
        mean_peak_distance = mean(peak_distances);
        std_peak_distance = std(peak_distances);

        PDmean_tmp(j) = mean_peak_distance;
        PDstd_tmp(j) = std_peak_distance;

%         clear T output omega1 y1;

    end

%     X(i, :, :) = XX;
    Y(i, :, :) = YY;
    PDmean(i, :) = PDmean_tmp;
    PDstd(i, :) = PDstd_tmp;

%     clear YY;

end

S = struct();
S.natural_period = natural_period;
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
% S.X = X;
S.PDmean = PDmean;
S.PDstd = PDstd;

date_string = datestr(clock());
filename = ['VanDerPol_ArnoldTongue_', date_string, '.mat'];
save(filename, '-struct', 'S');
