volume = inf;
% volume = 1e2;
% volume = 5e1;
% volume = 1e1;
% volume = 1e0;
% volume = 1e-1;
disp(volume);

omega = volume;

Ntrials = 1;

t0 = 0
tf = 100;
dt = 1e-2;

input_amplitude = 0.5;
input_period = 5;

Parameters();

% y0 = zeros(2,1);
y0 = [1; 1];

additive_forcing_func = @(t, x) AdditiveForcing(t, x, input_period, input_amplitude);
multiplicative_forcing_func = @(t, x) 0;

SDE = VanDerPol_Model(y0, omega, additive_forcing_func, multiplicative_forcing_func)
Nsteps = ceil((tf - t0) / dt);

[Paths, Times, Z] = SDE.simByEuler(Nsteps, 'DeltaTime', dt, 'NTRIALS', Ntrials);

T = Times;

output = Paths(:, 1);

figure();
plot(T, output);
title(['y(1): dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);


addpath('../');

% offset = find(T >= 100 * 60, 1);
% T = T(offset:end);
% output = output(offset:end);

[omega, y] = compute_fft(output, dt);
min_omega = 2 * pi * 0.01;
max_omega = 2 * pi * 1.0;
i1 = find(omega < min_omega, 1, 'last');
i2 = find(omega > max_omega, 1, 'first');
omega = omega(i1:i2);
y = y(i1:i2);

figure();
plot(omega ./ (2 * pi), abs(y) .^ 2);
title(['y(1) fft: dt=', num2str(dt), ' volume=', num2str(volume), ' amplitude=', num2str(input_amplitude), ' period=', num2str(input_period)]);
xlabel('frequency f');
ylabel('power |y|^2');

% filename = ['output/simulation_Ntrials=', int2str(Ntrials), ' dt=', num2str(dt), ' volume=', num2str(volume), ' offset=', num2str(TNF_offset), ' amplitude=', num2str(TNF_amplitude), ' period=', num2str(TNF_period), '.mat'];
% save(filename);
