% volume in liters
% time in minutes
function [T, P, Z, omega] = NFkB_Langevin(x0, tf, dt, volume, tnf_function, Ntrials, doPlot)

    if nargin < 6
        Ntrials = 1;
    end

    species_str = {'nNFkB', 'IkBm', 'cIkB', 'aIKK', 'iIKK'};
    Na = 6.022 * 1e23;
    omega = 1e-6 * Na * volume;
    SDE = NFkB_Langevin_Model(x0, omega, tnf_function);

    % initial time t0
%     t0 = 0;

    Nsteps = ceil(tf / dt);

    [P, T, Z] = SDE.simByEuler(Nsteps, 'DeltaTime', dt, 'NTRIALS', Ntrials);

    %Forcing = TNF_Forcing(T, omega, C_TNF, M_TNF);

    if (nargin >= 7 && doPlot) || (nargin < 7 && Ntrials == 1)
        X = P(:, :, 1);
        % plot trajectory
        figure();
        co = get(gca, 'ColorOrder');
        subplot(2, 1, 1);
        cla;
        hold on;
        for i = 1:size(X, 2)
            plot(T, X(:,i), 'Color', co(i, :));
        end
        legend(species_str);
        title('All trajectories');
        hold off;
        subplot(2, 1, 2);
        cla;
        hold on;
        plot(T, X(:,1));
        TNF = zeros(size(T));
        for i=1:length(T)
            TNF(i) = tnf_function(T(i), X(i, :));
        end
        yLimits = ylim();
        yMin = yLimits(1);
        yMax = yLimits(2);
        TNF = TNF - min(TNF);
        TNF = TNF / (max(TNF) - min(TNF));
        TNF = 0.5 * TNF + 0.25;
        TNF = (yMax - yMin) * TNF + yMin;
        plot(T, TNF, '--', 'Color', [0.5, 0.5, 0.5]);
        title('NFkB trajectory');
        hold off;
    end

end
