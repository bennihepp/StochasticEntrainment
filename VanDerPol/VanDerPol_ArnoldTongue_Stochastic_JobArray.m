% brutus
% bsub -n 1 -R "rusage[mem=1024]" -W 1:00 -J "job[1-624]" -o logs/VanDerPol_ArnoldTongue_Stochastic_JobArray_%I.out bash VanDerPol_ArnoldTongue_Stochastic_JobArray.sh "\$LSB_JOBINDEX"
%
% INDEX=14; bsub -n 8 -R "rusage[mem=4096]" -W 8:00 -R "span[ptile=8]" -o logs3/NFkB_Langevin_Arnold_Tongues_Brutus3_$INDEX.out bash NFkB_Langevin_Arnold_Tongues_Brutus3.sh $INDEX 8
%
% bsub -n 4 -w "done(56236368)" -R "rusage[mem=4096]" -W 8:00 -R "span[ptile=4]" -o logs3/process_arnold_tongues.out bash process_arnold_tongues3.sh 4

function VanDerPol_ArnoldTongue_Stochastic_JobArray(L)

    output_prefix = 'output/VanDerPol_ArnoldTongue_Stochastic_JobArray';

    addpath('../');

%     natural_period = 1/0.1065;

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


    l = L;
    display(['l=', int2str(l), ' out of ', int2str(num_of_iterations)]);
    i = per_ind(l);
    j = amp_ind(l);
    input_period = input_periods(i);
    input_amplitude = input_amplitudes(j);
    display(['input_period=', num2str(input_period), ' input_amplitude=', num2str(input_amplitude)]);

    additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
    multiplicative_forcing_func = @(t, x) 0;

    [TT, output] = VanDerPol_Run(Ntrials, t0, tf, dt, volume, additive_forcing_func, multiplicative_forcing_func);

%     X = output;

    offset_time = (tf - t0) / 5;
    offset_time = min(offset_time, 1000);
    offset = find(TT >= offset_time, 1);
    TT = TT(offset:end);
    output = output(offset:end, :);

    %% Fourier spectrum analysis
    omega = zeros(Ntrials, length(Omega));
    y = zeros(Ntrials, length(Omega));
    for m=1:Ntrials
        [omega1, y1] = compute_normalized_fft_truncated(output(:,m)', dt, 2*pi*min_frequency, 2*pi*max_frequency);
        omega(m, :) = omega1;
        y(m, :) = y1;
    end

%     Omega = omega;
    Y = mean(y, 1);
    Yabs = mean(abs(y), 1);

    %% Autocorrelation analysis
    corr = xcorr(mean(output,2) - mean(output(:)), 'unbiased');
    [~, locs] = findpeaks(corr);
    peak_distances = locs(2:end) - locs(1:end-1);
    mean_peak_distance = mean(peak_distances);
    std_peak_distance = std(peak_distances);

    PDmean = mean_peak_distance;
    PDstd = std_peak_distance;

    S = struct();
    S.Y = Y;
    S.Yabs = Yabs;
    % S.X = X;
    S.PDmean = PDmean;
    S.PDstd = PDstd; %#ok<STRNU>

    filename = [output_prefix, '_l=1', int2str(l), '.mat'];
    save(filename, '-struct', 'S');

end
