% S = load('output/NFkB_TNF_ArnoldTongue_BinarySearch_01-Sep-2014 14:24:46');
% S = load('output/NFkB_TNF_ArnoldTongue_BinarySearch_volume=Inf_population=0_09-Sep-2014 15:31:04');
% S = load('output_200/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=5e-13_population=0_10-Sep-2014 09:25:45');
% S = load('output_population_200/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=5e-13_population=1_10-Sep-2014 09:26:08');
% S = load('output_500/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=0_10-Sep-2014 23:16:20');
% S = load('output_population_500/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=1_11-Sep-2014 08:23:55');
% S = load('output_population_1000/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=1_13-Sep-2014 10:42:51');
S = load('output/NFkB_TNF_ArnoldTongue_BinarySearch_volume=Inf_population=0_19-Sep-2014 15:27:20');
% S = load('output_1000/NFkB_TNF_ArnoldTongue_BinarySearch_JobArray_volume=2e-11_population=0_19-Sep-2014 14:41:01');

input_periods = S.input_periods;
input_period_indices = 1:length(input_periods);
% i1 = find(S.input_periods >= 4.0, 1, 'first');
% i2 = find(S.input_periods <= 20.0, 1, 'last');
% input_period_indices = i1:i2;
% input_periods = S.input_periods(input_period_indices);

% input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:S.max_input_amplitude;
input_amplitudes = S.min_input_amplitude:S.input_amplitude_tolerance:S.max_input_amplitude;

Q = zeros(length(input_amplitudes), length(input_period_indices));
levels = 1;

for n=1:length(input_period_indices)
    i = input_period_indices(n);
    border = S.arnold_tongue_borders(i);
    j = find(input_amplitudes >= border, 1, 'first');
    Q(j:end, n) = 1;
end

if isfield(S, 'score_std')
    levels = 2;
    for n=1:length(input_period_indices)
        i = input_period_indices(n);
        border = S.arnold_tongue_borders(i);
        if isinf(border)
            continue;
        end
        stddev = S.score_std(i);
        if ~isinf(stddev) && ~isnan(stddev)
            j = find(input_amplitudes >= border, 1, 'first');
            Q(j:end, n) = 0.5;
        end
        j = find(input_amplitudes - stddev >= border, 1, 'first');
        Q(j:end, n) = 1.0;
    end
    levels = 2;
end

figure();
contourf(input_periods, input_amplitudes, Q, levels);
% title(['population arnold tongue for volume=', num2str(S.volume)]);
if isfield(S, 'population_average') && S.population_average
    title(['population arnold tongue for volume=', num2str(S.volume)]);
else
    title(['arnold tongue for volume=', num2str(S.volume)]);
end
xlabel('input period');
ylabel('input amplitude');
colorbar();
