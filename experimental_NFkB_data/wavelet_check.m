addpath([getenv('HOME'), '/Documents/MATLAB/wave_matlab']);

t = 1:0.01:(2*2*pi);
x = sin(t);
dt = mean(t(2:end) - t(1:end-1));
figure();
plot(t, x);

pad = 0;
dj = 0.1;
[wave,period,scale,coi] = wavelet(x, dt, pad, dj);
power = abs(wave) .^ 2;
phase = angle(wave);

figure();
contourf(t, period, power);
colorbar();
title('power');

% figure();
% contourf(t(2:end), period, phase(:, 2:end) - phase(:, 1:end-1));
% colorbar();
% title('phase');
