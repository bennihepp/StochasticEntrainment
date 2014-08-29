% input_periods = S.input_periods;
% input_period_indices = 1:length(input_periods);
i1 = find(S.input_periods >= 4.0, 1, 'first');
i2 = find(S.input_periods <= 20.0, 1, 'last');
input_period_indices = i1:i2;
input_periods = S.input_periods(input_period_indices);

% input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:S.max_input_amplitude;
input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:1.0;

Q = zeros(length(input_amplitudes), length(input_period_indices));
levels = 1;

for n=1:length(input_period_indices)
    i = input_period_indices(n);
    border = S.arnold_tongue_borders(i);
    j = find(input_amplitudes >= border, 1, 'first');
    Q(j:end, n) = 1;
end

if isfield(S, 'score_std')    
    for n=1:length(input_period_indices)
        i = input_period_indices(n);
        border = S.arnold_tongue_borders(i);
        if isinf(border)
            continue;
        end
        std = S.score_std(i);
        if ~isinf(std)
            j = find(input_amplitudes >= border, 1, 'first');
            Q(j:end, n) = 0.5;
        end
        j = find(input_amplitudes - std >= border, 1, 'first');
        Q(j:end, n) = 1.0;
    end
    levels = 2;
end

figure();
%contourf(input_periods, input_amplitudes, Q, levels);
surf(input_periods, input_amplitudes, double(Q),  'LineStyle', 'none');
xlim([min(input_periods), max(input_periods)]);
ylim([min(input_amplitudes), max(input_amplitudes)]);
view(0, 90);
% title(['population arnold tongue for volume=', num2str(S.volume)]);
title(['arnold tongue for volume=', num2str(S.volume)]);
xlabel('input period');
ylabel('input amplitude');
colorbar();
