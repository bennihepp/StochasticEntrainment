HEATMAP_TYPE = 'surface';

addpath('../');

S1 = load('output/NFkB_TNF_ArnoldTongue_BinarySearch_volume=Inf_population=0_09-Sep-2014 15:31:04');
% S2 = load('output_population_200/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=5e-13_population=1_10-Sep-2014 09:26:08');
% S2 = load('output_500/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=0_10-Sep-2014 23:16:20');
% S2 = load('output_1000/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=0_12-Sep-2014 22:26:37');
S2 = load('output_population_1000/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=1_13-Sep-2014 10:42:51');
% S2 = load('output_population_500/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=1_11-Sep-2014 08:23:55');

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

display(['area ratio: ', num2str(sum(Q2(:)) / sum(Q1(:)))]);
