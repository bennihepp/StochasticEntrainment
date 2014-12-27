HEATMAP_TYPE = 'matrix';
addpath('../');
addpath([getenv('HOME'), '/Documents/MATLAB/plotting']);

export_eps = true;
%export_eps = false;

% S1 = load('output/VanDerPol_ArnoldTongue_BinarySearch_01-Sep-2014 10:19:54');
% S2 = load('output/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_04-Sep-2014 23:58:25');
% S2 = load('output/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_30-Aug-2014 12:26:58');
% S2 = load('output/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population_average=true_30-Aug-2014 12:26:10');
% S1 = load('output/VanDerPol_ArnoldTongue_BinarySearch_volume=Inf_population=0_12-Sep-2014 09:11:13');
S1 = load('output_0.75/VanDerPol_ArnoldTongue_BinarySearch_volume=Inf_population=0_23-Oct-2014 12:54:48');
% S2 = load('output_200/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_11-Sep-2014 22:33:01');
% S2 = load('output_500/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_12-Sep-2014 12:04:41');
S2 = load('output_500_0.75/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_20-Sep-2014 12:17:56');
% S2 = load('output_population_500/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=1_12-Sep-2014 08:59:17');
% S2 = load('output_population_500_0.75/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=1_21-Sep-2014 10:55:20');

assert(all(S1.input_periods == S2.input_periods));

input_periods = S1.input_periods;
input_period_indices = 1:length(input_periods);

% input_amplitudes = S1.min_input_amplitude:S1.input_amplitude_tolerance:S1.max_input_amplitude;
input_amplitudes = S1.min_input_amplitude:S1.input_amplitude_tolerance:1.0;

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
    if ~isinf(S2.score_std(i)) && ~isnan(S2.score_std(i))
        % FIXME
        j = find(input_amplitudes >= border, 1, 'first') + 1;
    else
        j = find(input_amplitudes >= border, 1, 'first');
    end
%     j = find(input_amplitudes >= border, 1, 'first');
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
        filename = [export_eps_prefix(), 'vanderpol_comparison_arnold_tongue_population_volume1=', num2str(S1.volume), '_volume2=', num2str(S2.volume)];
    else
        filename = [export_eps_prefix(), 'vanderpol_comparison_arnold_tongue_volume1=', num2str(S1.volume), '_volume2=', num2str(S2.volume)];
    end
    save_plot(filename, h, width, height);
    %saveas(gcf(), filename, 'psc2');
end
