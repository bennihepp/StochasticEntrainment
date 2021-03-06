S = load('output/NFkB_TNF_ArnoldTongue_01-Sep-2014 13:31:07');

ENTRAINMENT_THRESHOLD = 0.9;

Q = S.scores' > ENTRAINMENT_THRESHOLD;

figure();
contourf(S.input_periods, S.input_amplitudes, Q, 1);
% title(['population arnold tongue for volume=', num2str(S.volume)]);
if isfield(S, 'population_average') && S.population_average
    title(['population arnold tongue for volume=', num2str(S.volume)]);
else
    title(['arnold tongue for volume=', num2str(S.volume)]);
end
xlabel('input period');
ylabel('input amplitude');
colorbar();
