input_periods = S.input_periods;
input_amplitudes = S.input_amplitudes;

Q = S.scores >= S.ENTRAINMENT_THRESHOLD;

Q(isinf(Q)) = max(Q(:));

figure();
contourf(input_periods, input_amplitudes, Q, 1);
title('arnold tongue');
xlabel('input period');
ylabel('input amplitude');
colorbar();
