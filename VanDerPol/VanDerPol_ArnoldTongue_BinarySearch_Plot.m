input_periods = S.input_periods;
% input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:S.max_input_amplitude;
input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:2.0;

Q = zeros(length(input_amplitudes), length(input_periods));

for i=1:length(input_periods)
    border = S.arnold_tongue_borders(i);
    j = find(input_amplitudes >= border, 1, 'first');
    Q(j:end, i) = 1;
end

figure();
contourf(input_periods, input_amplitudes, Q, 1);
title('arnold tongue');
xlabel('input period');
ylabel('input amplitude');
colorbar();
