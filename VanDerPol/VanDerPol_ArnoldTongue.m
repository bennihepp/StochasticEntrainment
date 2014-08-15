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

% input_amplitudes = 0.0:0.1:1.0;
% input_periods = 1:1:20;

min_frequency = 0.01;
max_frequency = 1.0;

q = t0:dt:tf;
[Omega, ~] = compute_fft_truncated(q, dt, 2*pi*min_frequency, 2*pi*max_frequency);
Y = zeros(length(input_periods), length(input_amplitudes), length(Omega));


parfor i=1:length(input_periods)
    display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
    input_period = input_periods(i);

    YY = zeros(length(input_amplitudes), length(Omega));

    for j=1:length(input_amplitudes)
        display(['j=', int2str(j), ' out of ', int2str(length(input_amplitudes))]);
        input_amplitude = input_amplitudes(j);

        additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
        multiplicative_forcing_func = @(t, x) 0;

        [T, output] = VanDerPol_Run(Ntrials, t0, tf, dt, omega, additive_forcing_func, multiplicative_forcing_func);

        offset_time = (tf - t0) / 5;
        offset_time = min(offset_time, 1000);
        offset = find(T >= offset_time, 1);
        T = T(offset:end);
        output = output(offset:end, :);

        [omega1, y1]= compute_fft_truncated(output, dt, 2*pi*min_frequency, 2*pi*max_frequency);
        min_omega = 2 * pi * min_frequency;
        max_omega = 2 * pi * max_frequency;
        i1 = find(omega1 < min_omega, 1, 'last');
        i2 = find(omega1 > max_omega, 1, 'first');
        omega1 = omega1(i1:i2);
        y1 = y1(i1:i2);

%         Omega = omega;
        YY(j, :) = y1;

%         clear T output omega1 y1;

    end

    Y(i, :, :) = YY;

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
S.Omega = Omega;
S.Y = Y;

filename = ['VanDerPol_ArnoldTongue_', date(), '.mat'];
save(filename, '-struct', 'S');
