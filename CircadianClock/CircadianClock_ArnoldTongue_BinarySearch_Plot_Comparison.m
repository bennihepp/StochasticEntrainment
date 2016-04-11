HEATMAP_TYPE = 'matrix';
addpath('../');
addpath([getenv('HOME'), '/Documents/MATLAB/plotting']);

export_eps = true;
% export_eps = false;

% S1 = load('output/CircadianClock_ArnoldTongue_BinarySearch_29-Aug-2014 11:28:54');
% S1 = load('output/CircadianClock_ArnoldTongue_BinarySearch_02-Sep-2014 16:30:40');
S1 = load('output/CircadianClock_ArnoldTongue_BinarySearch_09-Sep-2014 12:45:00');
% S1 = load('output_0.75/CircadianClock_ArnoldTongue_BinarySearch_23-Oct-2014 17:11:15');
% S2 = load('output/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=0_01-Sep-2014 16:33:13');
% S2 = load('output/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=1_01-Sep-2014 06:58:41');
% S2 = load('output_200/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=0_09-Sep-2014 23:14:44');
% S2 = load('output_1000/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=0_19-Sep-2014 14:37:29');
% S2 = load('output_1000_0.75/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=0_14-Oct-2014 14:52:36');
% S2 = load('output_population_200/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=1_09-Sep-2014 23:15:58');
S2 = load('output_population_1000/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=1_20-Oct-2014 17:11:10');
% S2 = load('output_population_1000_0.75/CircadianClock_ArnoldTongue_BinarySearch_JobArray_volume=1e-20_population=1_16-Oct-2014 10:38:08');


assert(all(S1.input_periods == S2.input_periods));

input_periods = S1.input_periods;
input_period_indices = 1:length(input_periods);

% i1 = find(S.input_periods >= 4.0, 1, 'first');
i1 = find(S1.input_periods >= 18, 1, 'first');
i2 = find(S1.input_periods <= 30.0, 1, 'last');
input_period_indices = i1:i2;
input_periods = S1.input_periods(input_period_indices);

% input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:S.max_input_amplitude;
input_amplitudes = S1.min_input_amplitude:S1.input_amplitude_tolerance:S1.max_input_amplitude;
input_amplitudes = S1.min_input_amplitude:S1.input_amplitude_tolerance:0.3;

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

if export_eps
    %set(gca(), 'FontSize', 20);
    width = 10;
    height = 6;
    fontSize = 0.5*0.5 * (width * height);
    h = prepare_plot(width, height, fontSize);
else
    figure();
end

arnold_tongue_colormap();
plot_heatmap(input_periods, input_amplitudes, Q, HEATMAP_TYPE, levels);
if ~export_eps
    title(['arnold tongue comparison for volume1=', num2str(S1.volume), ' and volume2=', num2str(S2.volume)]);
end
xlabel('input period');
ylabel('input amplitude');
% colorbar(...
%     'Location', 'SouthOutside', ...
%     'XTick', [0, 1, 2, 3], ...
%     'XTickLabel', { ...
%         'No entrainment', ...
%         ['Entrainment for volume=', num2str(S1.volume)], ...
%         ['Entrainment for volume=', num2str(S2.volume)], ...
%         'Entrainment for both', ...
%     } ...
% );

display(['area ratio: ', num2str(sum(Q2(:)) / sum(Q1(:)))]);

if export_eps
    if S2.population_average
        filename = [export_eps_prefix(), 'circadian_comparison_arnold_tongue_population_volume1=', num2str(S1.volume), '_volume2=', num2str(S2.volume)];
    else
        filename = [export_eps_prefix(), 'circadian_comparison_arnold_tongue_volume1=', num2str(S1.volume), '_volume2=', num2str(S2.volume)];
    end
    save_plot(filename, h, width, height);
%     saveas(gcf(), filename, 'psc2');
end
