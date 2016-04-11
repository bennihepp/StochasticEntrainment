initial_phases = linspace(0, 2*pi,10);
output_phases = zeros(size(initial_phases));
for i=1:length(initial_phases)
    initial_phase = initial_phases(i);
    VanDerPol_Batch;

    output = output - repmat(mean(output, 1), [size(output, 1), 1]);

    %% compute spectras
    addpath('../');
    [omega, ff_spectrum] = compute_normalized_fft_truncated(output, recordStep, 2*pi*min_frequency, 2*pi*max_frequency);

    [~, ind] = min(abs(omega ./ (2 * pi) - 1 ./ natural_period));
    output_phases(i) = angle(ff_spectrum(ind));

%     figure();
%     plot(output(1:1000));
end

% [~, ind] = min(output_phases);
% X = [initial_phases(ind:end), initial_phases(1:ind-1)];
% Y = [output_phases(ind:end), output_phases(1:ind-1)];

% figure();
% plot(initial_phases, output_phases);

% figure();
% hold on;
% plot(initial_phases(1:ind-1), output_phases(1:ind-1));
% plot(initial_phases(ind:end), output_phases(ind:end));

% figure();
% plot(X, Y);

X = initial_phases;
Y = output_phases;
Y(Y < Y(1)) = Y(Y < Y(1)) + 2*pi;
Y(end) = Y(end) + 2*pi;
X = X / (pi);
Y = (Y - Y(1)) / (pi);

xaxis_offset = 0.1;
yaxis_offset = 0.08;
width = 6;
height = 4;
fontSize = 0.5 * (width * height);
h = prepare_plot(width, height, fontSize);
hold on;
plot(X, Y, 'b', 'LineWidth', 2.0);
hold off;
hx = xlabel('input phase \phi');
hy = ylabel('output phase \theta');
set(gca, 'xlim', [0, 2], 'ylim', [0, 2]);
pos = get(gca, 'Position');
set(gca, 'Position', pos + [yaxis_offset, 0, -yaxis_offset, 0]);
% set(gca,...
%  'xlim',[0, 2],...
%  'xtick',[0, 0.5, 1, 1.5, 2.0],...
%  'xticklabel',{'0' 'p/2' 'p' '3/2p' '2p'},...
%  'ylim',[0, 2],...
%  'ytick',[0, 1, 2],...
%  'yticklabel',{'q0+0' 'q0+p' 'q0+2p'},...
%  'fontname','symbol',...
%  'fontsize',20);

format_ticks(gca, ...
    {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'}, ...
    {'\theta_{0}+0', '\theta_{0}+\pi', '\theta_{0}+2\pi'}, ...
    [0, 0.5, 1, 1.5, 2], ...
    [0, 1, 2]...
);

set(hx, 'Units', 'Normalized');
pos = get(hx, 'Position');
set(hx, 'Position', pos + [0, -xaxis_offset, 0]);

set(hy, 'Units', 'Normalized');
pos = get(hy, 'Position');
set(hy, 'Position', pos + [-yaxis_offset-0.05, 0, 0]);

save_plot([export_eps_prefix(), 'vanderpol_forcing_phase_synchronization'], h, width, height);
