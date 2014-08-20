addpath('../');

natural_period = 1/0.1065;

volume = 5e3;

if volume == inf
    Ntrials = 1;
    dt = 1e-1;
else
    Ntrials = 50;
    dt = 1e-1;
end

disp(['volume=', num2str(volume), ' Ntrials=', int2str(Ntrials), ' dt=', num2str(dt)]);

t0 = 0;
tf = 10000;


input_amplitudes = 0.0:0.1:2.0;
input_periods = 1:1.0:20;

% input_amplitudes = 0.0:0.2:1.5;
% input_periods = 2:2.0:30;

% input_amplitudes = 0.0:0.1:1.0;
% input_periods = 2:2:20;

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

% X = zeros(length(input_periods), length(input_amplitudes), length(T));
Y = zeros(length(input_periods), length(input_amplitudes), length(Omega));
Yabs = zeros(length(input_periods), length(input_amplitudes), length(Omega));
PDmean = zeros(length(input_periods), length(input_amplitudes));
PDstd = zeros(length(input_periods), length(input_amplitudes));


%for i=1:length(input_periods)
parfor i=1:length(input_periods)
    display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
    input_period = input_periods(i);

%     XX = zeros(length(input_amplitudes), length(T));
    YY = zeros(length(input_amplitudes), length(Omega));
    YYabs = zeros(length(input_amplitudes), length(Omega));
    PDmean_tmp = zeros(length(input_amplitudes), 1);
    PDstd_tmp = zeros(length(input_amplitudes), 1);

    for j=1:length(input_amplitudes)
        display(['i=', int2str(i), ' of ', int2str(length(input_periods)), ', j=', int2str(j), ' of ', int2str(length(input_amplitudes))]);
        input_amplitude = input_amplitudes(j);

        additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
        multiplicative_forcing_func = @(t, x) 0;

        [TT, output] = VanDerPol_Run(Ntrials, t0, tf, dt, volume, additive_forcing_func, multiplicative_forcing_func);

%         XX(j, :) = output;

	%size(output)
        offset_time = (tf - t0) / 5;
        offset_time = min(offset_time, 1000);
        offset = find(TT >= offset_time, 1);
        TT = TT(offset:end);
        output = output(offset:end, :);

        %% Fourier spectrum analysis
	%disp('q');
        omega = zeros(Ntrials, length(Omega));
        y = zeros(Ntrials, length(Omega));
	%size(omega)
	%size(y)
        for l=1:Ntrials
            [omega1, y1] = compute_normalized_fft_truncated(output(:,l)', dt, 2*pi*min_frequency, 2*pi*max_frequency);
		%disp('a');
		%size(omega1)
		%size(y1)
            omega(l, :) = omega1;
            y(l, :) = y1;
        end
%         mean_y = mean(y, 1);
%         mean_omega = mean(omega, 1);

%         Omega = omega;
		%disp('b');
	%size(mean(y,1))
	%size(mean(abs(y), 1))
        YY(j, :) = mean(y, 1);
        YYabs(j, :) = mean(abs(y), 1);

        %% Autocorrelation analysis
		%disp('c');
	%size(mean(output,2))
	%size(mean(output(:)))
        corr = xcorr(mean(output,2) - mean(output(:)), 'unbiased');
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
    Yabs(i, :, :) = YYabs;
    PDmean(i, :) = PDmean_tmp;
    PDstd(i, :) = PDstd_tmp;

%     clear YY;

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
