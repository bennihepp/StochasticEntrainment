HEATMAP_TYPE = 'surface';

addpath('../');

S1 = load('output/CircadianClock_ArnoldTongue_BinarySearch_29-Aug-2014 11:28:54');
% S1 = load('output/CircadianClock_ArnoldTongue_BinarySearch_02-Sep-2014 16:30:40');
% S2 = load('output/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=0_01-Sep-2014 16:33:13');
S2 = load('output_200/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=0_09-Sep-2014 09:53:49');
% S2 = load('output/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=1_01-Sep-2014 06:58:41');
% S2 = load('output_population_200/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=1_09-Sep-2014 08:15:38');

assert(all(S1.input_periods == S2.input_periods));

input_periods = S1.input_periods;
input_period_indices = 1:length(input_periods);
% i1 = find(S.input_periods >= 4.0, 1, 'first');
% i2 = find(S.input_periods <= 20.0, 1, 'last');
% input_period_indices = i1:i2;
% input_periods = S.input_periods(input_period_indices);

% input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:S.max_input_amplitude;
input_amplitudes = S1.min_input_amplitude:S1.input_amplitude_tolerance:S1.max_input_amplitude;

Q1 = zeros(length(input_amplitudes), length(input_period_indices));
Q2 = zeros(length(input_amplitudes), length(input_period_indices));
levels = 3;

for n=1:length(input_period_indices)
    i = input_period_indices(n);
    % S1
    border = S1.arnold_tongue_borders(i);
    j = find(input_amplitudes >= border, 1, 'first');
    Q1(j:end, n) = 1;
    % S2
    border = S2.arnold_tongue_borders(i);
%     if ~isinf(S2.score_std(i))
%         j = find(input_amplitudes >= border, 1, 'first') + 1;
%     else
%         j = find(input_amplitudes >= border, 1, 'first');
%     end
    j = find(input_amplitudes >= border, 1, 'first');
    Q2(j:end, n) = 1;
end

Q = Q1 + 2 * Q2;

figure();
plot_heatmap(input_periods, input_amplitudes, Q, HEATMAP_TYPE, levels);
title(['arnold tongue comparison for volume1=', num2str(S1.volume), ' and volume2=', num2str(S2.volume)]);
xlabel('input period');
ylabel('input amplitude');
colorbar(...
    'Location', 'SouthOutside', ...
    'XTick', [0, 1, 2, 3], ...
    'XTickLabel', { ...
        'No entrainment', ...
        ['Entrainment for volume=', num2str(S1.volume)], ...
        ['Entrainment for volume=', num2str(S2.volume)], ...
        'Entrainment for both', ...
    } ...
);
