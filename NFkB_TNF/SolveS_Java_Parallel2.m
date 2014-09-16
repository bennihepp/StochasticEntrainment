% volume in liters
% time in minutes
function [T, P, omega] = SolveS_Java_Parallel2(x0, tf, dt, volume, ...
                            input_offset, input_amplitude, input_frequency, Ntrials, ...
                            recordStep, seed, printProgress)

    DEFAULT_NUM_OF_THREADS = 8;

    NthreadsStr = getenv('NUM_OF_THREADS');
    if ~isempty(NthreadsStr)
        Nthreads = str2double(NthreadsStr);
    else
        Nthreads = DEFAULT_NUM_OF_THREADS;
    end

    if nargin < 9
        recordStep = dt;
    end
    if nargin < 8
        Ntrials = 1;
    end

    if nargin < 10
        seed = randi([0, 2.^31-2]);
    end

    if nargin < 11
        printProgress = false;
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
    sde = ch.ethz.bhepp.sdesolver.models.NFkBSpikySde(volume, inputFunction);

%     initial time t0
    t0 = 0;

    if size(x0, 2) == 1
        x0 = repmat(x0, 1, Ntrials);
    end

    javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
    javaaddpath(JavaLangevinModel_path);
    try
        ch.ethz.bhepp.sdesolver.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display(e);
    end

%     for n=1:Ntrials
    javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
    javaaddpath(JavaLangevinModel_path);
    try
        ch.ethz.bhepp.sdesolver.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display(e);
    end

%     omega = sde.getOmega();
    sdeFactory = ch.ethz.bhepp.sdesolver.models.NFkBSpikySdeFactory(volume, inputFunction);
    stepperFactory = ch.ethz.bhepp.sdesolver.EulerMaruyamaStepperFactory(sdeFactory, dt, seed);
    positiveStateHook = ch.ethz.bhepp.sdesolver.EnsurePositiveStateHook();
    stepperFactory.addStepHook(positiveStateHook);
%     boundedStateHook = ch.ethz.bhepp.sdesolver.EnsureBoundedStateHook(0, 0.0, 1.0);
%     stepperFactory.addStepHook(boundedStateHook);
    parallelSolver = ch.ethz.bhepp.sdesolver.ParallelSdeSolver(stepperFactory, Nthreads);

    x0Matrix = ch.ethz.bhepp.sdesolver.MatrixHelper.createZeroDoubleMatrix2D(size(x0, 1), Ntrials);
    for i=1:size(x0, 1)
        for n=1:Ntrials
            x0Matrix.set(i-1, n-1, x0(i, n));
        end
    end
    sdeSolutions = parallelSolver.solve(Ntrials, 0, x0Matrix, tf, recordStep, printProgress);

    for n=1:Ntrials
        sdeSolution = sdeSolutions.get(n - 1);
        TArray = sdeSolution.getTArray();
        XArray = sdeSolution.getXArray();
        if n == 1
            numOfTimeSteps = length(TArray);
            T = zeros(Ntrials, numOfTimeSteps);
            P = zeros(Ntrials, numOfTimeSteps, sde.getDriftDimension());
        end
        T(n, :) = TArray;
        P(n, :, :) = XArray;
    end
    T = T(1, :);

    inputFunction = ch.ethz.bhepp.sdesolver.SinusoidalFunction(input_offset, input_amplitude, input_frequency);
    sde = ch.ethz.bhepp.sdesolver.models.NFkBSpikySde(volume, inputFunction);
    omega = sde.getOmega();
end
