% S = load('output_0.9/ArnoldTongue_BinarySearch_28-Dec-2014 15:30:31');
% S = load('output_0.75/ArnoldTongue_BinarySearch_02-Jan-2015 11:43:02');
% S = load('output_0.75/ArnoldTongue_BinarySearch_02-Jan-2015 16:02:56');
% S = load('output_population_100/ArnoldTongue_BinarySearch_JobArray_population=1_01-Jan-2015 18:04:21');
% S = load('output_100_0.75/ArnoldTongue_BinarySearch_JobArray_population=0_04-Jan-2015 21:31:47');
S = load('output_population_100_0.75/ArnoldTongue_BinarySearch_JobArray_population=1_05-Jan-2015 10:06:47');

input_periods = S.input_periods;
input_period_indices = 1:length(input_periods);
% i1 = find(S.input_periods >= 4.0, 1, 'first');
% i2 = find(S.input_periods <= 20.0, 1, 'last');
% i1 = 1;
% i2 = find(S.input_periods <= 35.0, 1, 'last');
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
end

figure();
contourf(input_periods, input_amplitudes, Q, levels);
% title(['population arnold tongue for volume=', num2str(S.volume)]);
if isfield(S, 'population_average') && S.population_average
    title(['population arnold tongue']);
else
    title(['arnold tongue']);
end
xlabel('input period');
ylabel('input amplitude');
colorbar();
