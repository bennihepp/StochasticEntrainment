% volume in liters
% time in minutes
function [T, P, omega] = SolveS_Java(x0, tf, dt, volume, ...
                            input_offset, input_amplitude, input_frequency, Ntrials, ...
                            doPlot, recordStep, seed)

    if nargin < 10
        recordStep = dt;
    end
    if nargin < 9
        doPlot = true;
    end
    if nargin < 8
        Ntrials = 1;
    end

    if nargin < 11
        seed = randi([0, 2.^31-2]);
%         seed = 4942634935159178;
        disp(['Java seed=', int2str(seed)]);
    end

    if Ntrials > 1 && volume == inf && size(x0, 1) == 1
        Ntrials = 1;
    end

%         global JavaLangevinModel
%         JavaLangevinModel_path = '/cluster/home02/misc/bhepp/StochasticEntrainment/JavaLangevinModel/bin';
%         import_java();
    javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
    if ~exist('JavaLangevinModel_path', 'var') || ~ischar(JavaLangevinModel_path) %#ok<NODEF>
        NFkB_Langevin_Model_path = '../JavaLangevinModel/bin';
    end
%         disp(NFkB_Langevin_Model_path);
    javaaddpath(NFkB_Langevin_Model_path);
    try
        ch.ethz.bhepp.sdesolver.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display(e);
    end

    species_str = {'M', 'P_0', 'P_1', 'P_2', 'P_N'};
    inputFunction = ch.ethz.bhepp.sdesolver.SinusoidalFunction(input_offset, input_amplitude, input_frequency);
    sde = ch.ethz.bhepp.sdesolver.models.CircadianClockDrosophilaSde(volume, inputFunction);
    omega = sde.getOmega();
    stepper = ch.ethz.bhepp.sdesolver.EulerMaruyamaStepper(sde, dt, seed);
    positiveStateHook = ch.ethz.bhepp.sdesolver.EnsurePositiveStateHook();
    stepper.addStepHook(positiveStateHook);
%     boundedStateHook = ch.ethz.bhepp.sdesolver.EnsureBoundedStateHook(0, 0.0, 1.0);
%     stepper.addStepHook(boundedStateHook);
    solver = ch.ethz.bhepp.sdesolver.SdeSolver(stepper);

    % initial time t0
    %     t0 = 0;

    import ch.ethz.bhepp.sdesolver.*;

    x0Matrix = MatrixHelper.createZeroDoubleMatrix1D(size(x0, 1));
    if size(x0, 2) == 1
        x0 = repmat(x0, 1, Ntrials);
    end

    for i=1:size(x0, 1)
        x0Matrix.set(i-1, x0(i, 1));
    end
    sdeSolution = solver.solve(0, x0Matrix, tf, recordStep);
    TArray = sdeSolution.getTArray();
    XArray = sdeSolution.getXArray();
    T = TArray;
    P = zeros(Ntrials, size(XArray, 1), size(XArray, 2));
    for n=1:Ntrials
        if n > 1
            for i=1:size(x0, 1)
                x0Matrix.set(i-1, x0(i, n));
            end
            sdeSolution = solver.solve(0, x0Matrix, tf, recordStep);
%                 TArray = sdeSolution.getTArray();
            XArray = sdeSolution.getXArray();
        end
        P(n, :, :) = XArray;
    end

    %Forcing = TNF_Forcing(T, omega, C_TNF, M_TNF);

    if (nargin >= 7 && doPlot) || (nargin < 7 && Ntrials == 1)
        X = squeeze(P(1, :, :));
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
            XMatrix = MatrixHelper.createDoubleMatrix1D(X(i, :));
            TNF(i) = inputFunction.evaluate(T(i), XMatrix);
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
