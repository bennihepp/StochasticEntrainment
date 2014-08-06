t = 1:0.1:(100*2*pi);
x = sin(t);
dt = mean(t(2:end) - t(1:end-1));
figure();
plot(t, x);

[omega1, y1] = compute_fft(x, dt);
freq1 = omega1 / (2*pi);

i1 = find(omega1 >= 0.1, 1);
i2 = find(omega1 <= 10, 1, 'last');
figure();
plot(omega1(i1:i2), abs(y1(i1:i2)) .^ 2);


dfreq = 0.05;
disp('average around frequency 1');
freq = 1;
omega = omega1;
y = y1;
mask = (omega >= freq-dfreq) & (omega <= freq+dfreq);
q1 = sum(abs(mean(y(mask))) .^ 2);
disp([' complex average:', num2str(q1)]);



t = 1:0.1:(10*100*2*pi);
x = sin(t);
dt = mean(t(2:end) - t(1:end-1));
figure();
plot(t, x);

[omega2, y2] = compute_fft(x, dt);
freq2 = omega2 / (2*pi);

i1 = find(omega2 >= 0.1, 1);
i2 = find(omega2 <= 10, 1, 'last');
figure();
plot(omega2(i1:i2), abs(y2(i1:i2)) .^ 2);


dfreq = 0.05;
disp('average around frequency 1');
freq = 1;
omega = omega2;
y = y2;
mask = (omega >= freq-dfreq) & (omega <= freq+dfreq);
q1 = sum(abs(mean(y(mask))) .^ 2);
disp([' complex average:', num2str(q1)]);


% dfreq = 0.05;
% disp('average around frequency 1');
% freq = 1;
% omega = omega1;
% y = y1;
% mask = (omega >= freq-dfreq) & (omega <= freq+dfreq);
% q1 = sum(abs(mean(y(mask))) .^ 2);
% disp([' complex average:', num2str(q1)]);
% if Ntrials > 1
%     w1 = sum(mean_abs_omega(mask));
%     disp([' absolute average:', num2str(w1)]);
% end
% disp('average around frequency 0.50');
% freq = 0.50;
% freq2 = freq;
% mask = (mean_freq >= freq-dfreq) & (mean_freq <= freq+dfreq);
% q2 = sum(abs_mean_omega(mask));
% disp([' complex average:', num2str(q2)]);
% if Ntrials > 1
%     w2 = sum(mean_abs_omega(mask));
%     disp([' absolute average:', num2str(w2)]);
% end
% disp(['ratios ', num2str(freq2), ' to ', num2str(freq1)]);
% disp([' complex: ', num2str(q2 / q1)]);
% if Ntrials > 1
%     disp([' absolute: ', num2str(w2 / w1)]);
% end
