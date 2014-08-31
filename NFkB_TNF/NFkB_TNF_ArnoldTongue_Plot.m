ENTRAINMENT_THRESHOLD = 0.9;

Q = S.scores' > ENTRAINMENT_THRESHOLD;

figure();
contourf(input_periods, input_amplitudes, Q, 1);
% title(['population arnold tongue for volume=', num2str(S.volume)]);
if isfield(S, 'population_average') && S.population_average
    title(['population arnold tongue for volume=', num2str(S.volume)]);
else
    title(['arnold tongue for volume=', num2str(S.volume)]);
end
xlabel('input period');
ylabel('input amplitude');
colorbar();
