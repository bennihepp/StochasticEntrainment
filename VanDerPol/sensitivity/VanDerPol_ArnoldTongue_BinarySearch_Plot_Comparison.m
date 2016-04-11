HEATMAP_TYPE = 'matrix';
addpath('../../');
addpath('../../plotting');

export_eps_prefix = @() 'figures/';
export_eps = true;
% export_eps = false;

% S1 = load('output_inf_thres_0.75_B_10_d_2/VanDerPol_ArnoldTongue_BinarySearch__B_10_d_2_volume=Inf_population=0_21-Apr-2015 23:16:55');
% S2 = load('output_500_thres_0.75_B_10_d_2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_21-Apr-2015 08:42:15');
% % S2 = load('output_population_500_thres_0.75_B_10_d_2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=1_21-Apr-2015 08:42:23');

% S1 = load('output_inf_thres_0.75_B_9_d_2/VanDerPol_ArnoldTongue_BinarySearch__B_9_d_2_volume=Inf_population=0_21-Apr-2015 23:13:06');
% S2 = load('output_500_thres_0.75_B_9_d_2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_21-Apr-2015 08:39:04');
% % S2 = load('output_population_500_thres_0.75_B_9_d_2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=1_21-Apr-2015 08:40:00');

% S1 = load('output_inf_thres_0.75_B_11_d_2/VanDerPol_ArnoldTongue_BinarySearch__B_11_d_2_volume=Inf_population=0_22-Apr-2015 10:30:54');
% S2 = load('output_500_thres_0.75_B_11_d_2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_20-Apr-2015 23:24:56');
% % S2 = load('output_population_500_thres_0.75_B_11_d_2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=1_21-Apr-2015 08:24:06');

% S1 = load('output_inf_thres_0.75_B_10_d_1.8/VanDerPol_ArnoldTongue_BinarySearch__B_10_d_1.8_volume=Inf_population=0_21-Apr-2015 23:21:45.mat');
% S2 = load('output_500_thres_0.75_B_10_d_1.8/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_21-Apr-2015 08:26:58');
% % S2 = load('output_population_500_thres_0.75_B_10_d_1.8/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=1_21-Apr-2015 15:01:24');

S1 = load('output_inf_thres_0.75_B_10_d_2.2/VanDerPol_ArnoldTongue_BinarySearch__B_10_d_2.2_volume=Inf_population=0_22-Apr-2015 10:30:42.mat');
S2 = load('output_500_thres_0.75_B_10_d_2.2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=0_21-Apr-2015 17:31:25');
% S2 = load('output_population_500_thres_0.75_B_10_d_2.2/VanDerPol_ArnoldTongue_BinarySearch_JobArray_volume=5000_population=1_21-Apr-2015 14:56:16');

assert(S1.parameterB == S2.parameterB);
assert(S1.parameterD == S2.parameterD);

parameterB = S2.parameterB;
parameterD = S2.parameterD;

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

%     if isfinite(S2.score_std(i))
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
        filename = [export_eps_prefix(), 'vanderpol_sensitivity_B_', num2str(parameterB), '_d_', num2str(parameterD), '_comparison_arnold_tongue_population_volume1=', num2str(S1.volume), '_volume2=', num2str(S2.volume)];
    else
        filename = [export_eps_prefix(), 'vanderpol_sensitivity_B_', num2str(parameterB), '_d_', num2str(parameterD), '_comparison_arnold_tongue_volume1=', num2str(S1.volume), '_volume2=', num2str(S2.volume)];
    end
%     display(filename);
    save_plot(filename, h, width, height);
    saveas(gcf(), [filename, '.fig']);
%     saveas(gcf(), filename, 'psc2');
end
