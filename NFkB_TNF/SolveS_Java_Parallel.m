% volume in liters
% time in minutes
function [T, P, omega] = SolveS_Java_Parallel(x0, tf, dt, volume, ...
                            input_offset, input_amplitude, input_frequency, Ntrials, ...
                            recordStep)

    if nargin < 9
        recordStep = dt;
    end
    if nargin < 8
        Ntrials = 1;
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
        ch.ethz.bhepp.sdesolver.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display(e);
    end

    inputFunction = ch.ethz.bhepp.sdesolver.SinusoidalFunction(input_offset, input_amplitude, input_frequency);
    sde = ch.ethz.bhepp.sdesolver.models.CircadianClockDrosophilaSde(volume, inputFunction);

%     initial time t0
    t0 = 0;
    numOfTimeSteps = length(t0:recordStep:tf);

    if size(x0, 2) == 1
        x0 = repmat(x0, 1, Ntrials);
    end

    T = zeros(Ntrials, numOfTimeSteps);
    P = zeros(Ntrials, numOfTimeSteps, sde.getDriftDimension);
    par_seeds = zeros(Ntrials, 1);
    for n=1:Ntrials
        par_seeds(n) = randi([0, 2.^31-2]);
    end

    javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
    javaaddpath(JavaLangevinModel_path);
    try
        ch.ethz.bhepp.sdesolver.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display(e);
    end

    parfor n=1:Ntrials
        javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
        javaaddpath(JavaLangevinModel_path);
        try
            ch.ethz.bhepp.sdesolver.SinusoidalFunction(0.1, 0.05, 1/120);
        catch e
            display(e);
        end

        inputFunction = ch.ethz.bhepp.sdesolver.SinusoidalFunction(input_offset, input_amplitude, input_frequency);
        sde = ch.ethz.bhepp.sdesolver.models.CircadianClockDrosophilaSde(volume, inputFunction);
%         omega = sde.getOmega();
        stepper = ch.ethz.bhepp.sdesolver.EulerMaruyamaStepper(sde, dt, par_seeds(n));
        positiveStateHook = ch.ethz.bhepp.sdesolver.EnsurePositiveStateHook();
        stepper.addStepHook(positiveStateHook);
%         boundedStateHook = ch.ethz.bhepp.sdesolver.EnsureBoundedStateHook(0, 0.0, 1.0);
%         stepper.addStepHook(boundedStateHook);
        solver = ch.ethz.bhepp.sdesolver.SdeSolver(stepper);

        x0Matrix = ch.ethz.bhepp.sdesolver.MatrixHelper.createZeroDoubleMatrix1D(size(x0, 1));
        for i=1:size(x0, 1)
            x0Matrix.set(i-1, x0(i, n));
        end
        sdeSolution = solver.solve(0, x0Matrix, tf, recordStep);
        TArray = sdeSolution.getTArray();
        XArray = sdeSolution.getXArray();

        T(n, :) = TArray;
        P(n, :, :) = XArray;
    end
    T = T(1, :);

    inputFunction = ch.ethz.bhepp.sdesolver.SinusoidalFunction(input_offset, input_amplitude, input_frequency);
    sde = ch.ethz.bhepp.sdesolver.models.CircadianClockDrosophilaSde(volume, inputFunction);
    omega = sde.getOmega();
end
