HEATMAP_TYPE = 'matrix';
addpath('../');
addpath([getenv('HOME'), '/Documents/MATLAB/plotting']);

% export_eps = true;
export_eps = false;

% S1 = load('output_0.9/ArnoldTongue_BinarySearch_28-Dec-2014 14:24:07');
% S2 = load('output_population_50/ArnoldTongue_BinarySearch_JobArray_population=1_01-Jan-2015 18:03:00');

% S2 = load('output_50/ArnoldTongue_BinarySearch_JobArray_population=0_28-Dec-2014 19:37:22');
% S2 = load('output_100/ArnoldTongue_BinarySearch_JobArray_population=0_04-Jan-2015 13:45:00');
% S2 = load('output_population_100/ArnoldTongue_BinarySearch_JobArray_population=1_29-Dec-2014 22:01:58');
% S2 = load('output_population_100/ArnoldTongue_BinarySearch_JobArray_population=1_01-Jan-2015 18:04:21');

% S1 = load('output_0.75/ArnoldTongue_BinarySearch_02-Jan-2015 16:02:56');
% S1 = load('output_0.9/ArnoldTongue_BinarySearch_28-Dec-2014 15:30:31');
% S2 = load('output_population_500_0.75/ArnoldTongue_BinarySearch_JobArray_population=1_06-Jan-2015 10:09:32');
% S2 = load('output_population_500/ArnoldTongue_BinarySearch_JobArray_population=1_16-Jan-2015 13:05:35');
S2 = load('output_500_0.75/ArnoldTongue_BinarySearch_JobArray_population=0_11-Jan-2015 12:40:59');
% S2 = load('output_500/ArnoldTongue_BinarySearch_JobArray_population=0_19-Jan-2015 11:42:18');

% S1 = load('output_100_0.75/ArnoldTongue_BinarySearch_JobArray_population=0_04-Jan-2015 21:31:47');
% S2 = load('output_population_100_0.75/ArnoldTongue_BinarySearch_JobArray_population=1_05-Jan-2015 10:06:47');


assert(all(S1.input_periods == S2.input_periods));

input_periods = S1.input_periods;
input_period_indices = 1:length(input_periods);

% i1 = find(S.input_periods >= 4.0, 1, 'first');
% i1 = 1;
% i2 = find(S1.input_periods <= 35.0, 1, 'last');
% input_period_indices = i1:i2;
% input_periods = S1.input_periods(input_period_indices);

% input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:S.max_input_amplitude;
input_amplitudes = S1.min_input_amplitude:S1.input_amplitude_tolerance:S1.max_input_amplitude;
% input_amplitudes = S1.min_input_amplitude:S1.input_amplitude_tolerance:0.6;

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
    title(['arnold tongue comparison']);
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
        filename = [export_eps_prefix(), 'leloup_goldbeter_circadian_comparison_arnold_tongue_population'];
    else
        filename = [export_eps_prefix(), 'leloup_goldbeter_circadian_comparison_arnold_tongue'];
    end
    save_plot(filename, h, width, height);
%     saveas(gcf(), filename, 'psc2');
end
