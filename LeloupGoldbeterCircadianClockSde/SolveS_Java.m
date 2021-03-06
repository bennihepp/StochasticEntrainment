% volume in liters
% time in minutes
function [T, P, omega, sde] = SolveS_Java(x0, tf, dt, volume, ...
                            input_offset, input_amplitude, input_frequency, initial_phase, Ntrials, ...
                            recordStep, seed)

    if nargin < 10
        recordStep = dt;
    end
    if nargin < 9
        Ntrials = 1;
    end

    if nargin < 11
        seed = randi([0, 2.^31-2]);
%         seed = 4942634935159178;
%         disp(['Java seed=', int2str(seed)]);
    end

    if Ntrials > 1 && volume == inf && size(x0, 1) == 1
        Ntrials = 1;
    end

%         global JavaLangevinModel
%         JavaLangevinModel_path = '/cluster/home02/misc/bhepp/StochasticEntrainment/JavaLangevinModel/bin';
%         import_java();
    javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
    if ~exist('JavaLangevinModel_path', 'var') || ~ischar(JavaLangevinModel_path) %#ok<NODEF>
        JavaLangevinModel_path = '../JavaLangevinModel/bin';
    end
%         disp(JavaLangevinModel_path);
    javaaddpath(JavaLangevinModel_path);
    try
        ch.ethz.bhepp.utils.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display('Java error');
        display(e);
    end

    inputFunction = ch.ethz.bhepp.utils.SinusoidalFunction(input_offset, input_amplitude, input_frequency, initial_phase);
    sde = ch.ethz.bhepp.sdesolver.models.LeloupGoldbeterCircadianClockSde(volume, inputFunction);

    omega = sde.getOmega();
    stepper = ch.ethz.bhepp.sdesolver.EulerMaruyamaStepper(sde, dt, seed);
    positiveStateHook = ch.ethz.bhepp.sdesolver.EnsurePositiveStateHook();
    stepper.addStepHook(positiveStateHook);
%     boundedStateHook = ch.ethz.bhepp.sdesolver.EnsureBoundedStateHook(0, 0.0, 1.0);
%     stepper.addStepHook(boundedStateHook);
    solver = ch.ethz.bhepp.sdesolver.SdeSolver(stepper);

%     initial time t0
    t0 = 0;
    numOfTimeSteps = length(t0:recordStep:tf);

    x0Matrix = ch.ethz.bhepp.utils.MatrixHelper.createZeroDoubleMatrix1D(size(x0, 1));
    if size(x0, 2) == 1
        x0 = repmat(x0, 1, Ntrials);
    end

    for i=1:size(x0, 1)
        x0Matrix.set(i-1, x0(i, 1));
    end
    solution = solver.solve(0, x0Matrix, tf, recordStep);
    TArray = solution.getTArray();
    XArray = solution.getXArray();
%     T = TArray;
%     P = zeros(Ntrials, size(XArray, 1), size(XArray, 2));
    T = zeros(Ntrials, numOfTimeSteps);
    P = zeros(Ntrials, numOfTimeSteps, sde.getDriftDimension());

    for n=1:Ntrials
        display(['  iteration ', int2str(n), ' out of ', int2str(Ntrials)]);
        if n > 1
            for i=1:size(x0, 1)
                x0Matrix.set(i-1, x0(i, n));
            end
            solution = solver.solve(0, x0Matrix, tf, recordStep);
            TArray = solution.getTArray();
            XArray = solution.getXArray();
        end

        T(n, :) = TArray;
        P(n, :, :) = XArray;
    end
    T = T(1, :);

end
