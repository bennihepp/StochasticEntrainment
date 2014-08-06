per75 = nfkb_data.per75;
dt = per75.dt;
traces = per75.traces;
DT = mean(dt(2:end) - dt(1:end-1));

traces30 = traces(:,1:30);
% mask = ones(1, 30);
% mask(15) = 0;
% mask(20) = 0;
% mask(21) = 0;
% mask(23) = 0;
% mask(24) = 0;
% mask = repmat(mask, size(traces30, 1), 1);
% 
% traces = traces30(mask == 1);

traces30(:, 15) = [];
traces30(:, 19) = [];
traces30(:, 19) = [];
traces30(:, 20) = [];
traces30(:, 20) = [];

tr = traces(:, 4);
[omega, y] = compute_fft(tr, DT);
i1 = find(freq > 1/200, 1);
i2 = find(freq < 1/30, 1, 'last');
figure();
plot(freq(i1:i2), abs(y(i1:i2)).^2)
