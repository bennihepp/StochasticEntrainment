input_periods = S.input_periods;
input_amplitudes = S.input_amplitudes;

Q = S.scores >= S.ENTRAINMENT_THRESHOLD;

Q(isinf(Q)) = max(Q(:));

Q = Q';
figure();
% contourf(input_periods, input_amplitudes, Q, 1);
surf(input_periods, input_amplitudes, double(Q),  'LineStyle', 'none');
xlim([min(input_periods), max(input_periods)]);
ylim([min(input_amplitudes), max(input_amplitudes)]);
view(0, 90);
title('arnold tongue');
xlabel('input period');
ylabel('input amplitude');
colorbar();
