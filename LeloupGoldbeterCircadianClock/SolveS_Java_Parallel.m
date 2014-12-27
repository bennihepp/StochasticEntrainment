% volume in liters
% time in minutes
function [T, P, model] = SolveS_Java_Parallel(x0, t0, tf, recordStep, omega, ...
                            input_offset, input_amplitude, ...
                            input_frequency, initial_phase, Ntrials, ...
                            printMessages, seed)

    DEFAULT_NUM_OF_THREADS = 8;

    NthreadsStr = getenv('NUM_OF_THREADS');
    if ~isempty(NthreadsStr)
        Nthreads = str2double(NthreadsStr);
    else
        Nthreads = DEFAULT_NUM_OF_THREADS;
    end

    if nargin < 10
        Ntrials = 1;
    end

    if nargin < 11
        printMessages = false;
    end

    if nargin < 12
        seed = randi([0, 2.^31-2]);
    end

%     if Ntrials > 1 && volume == inf && size(x0, 1) == 1
%         Ntrials = 1;
%     end

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
        display(e);
    end

    inputFunction = ch.ethz.bhepp.utils.SinusoidalFunction(input_offset, input_amplitude, input_frequency, initial_phase);
    model = ch.ethz.bhepp.ssasolver.models.LeloupGoldbeterCircadianClock(omega, inputFunction);

    if size(x0, 2) == 1
        x0 = repmat(x0, 1, Ntrials);
    end

    javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
    javaaddpath(JavaLangevinModel_path);
    try
        ch.ethz.bhepp.utils.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display(e);
    end

%     for n=1:Ntrials
    javaaddpath([getenv('HOME'), '/local/lib/java/colt.jar']);
    javaaddpath(JavaLangevinModel_path);
    try
        ch.ethz.bhepp.utils.SinusoidalFunction(0.1, 0.05, 1/120);
    catch e
        display(e);
    end

%     omega = sde.getOmega();
    modelFactory = ch.ethz.bhepp.ssasolver.models.LeloupGoldbeterCircadianClockFactory(omega, inputFunction);
    stepperFactory = ch.ethz.bhepp.ssasolver.SimpleSSAStepperFactory(modelFactory, seed);

    stepperFactory = ch.ethz.bhepp.ssasolver.SimpleSSAStepperFactory(modelFactory, seed);
    parallelSolver = ch.ethz.bhepp.ssasolver.ParallelSSASolver(stepperFactory, Nthreads);

    x0Matrix = ch.ethz.bhepp.utils.MatrixHelper.createZeroDoubleMatrix2D(size(x0, 1), Ntrials);
    for i=1:size(x0, 1)
        for n=1:Ntrials
            x0Matrix.set(i-1, n-1, x0(i, n));
        end
    end
    sdeSolutions = parallelSolver.solve(Ntrials, t0, x0Matrix, tf, recordStep, printMessages);

    for n=1:Ntrials
        sdeSolution = sdeSolutions.get(n - 1);
        TArray = sdeSolution.getTArray();
        XArray = sdeSolution.getXArray();
        if n == 1
            numOfTimeSteps = length(TArray);
            T = zeros(Ntrials, numOfTimeSteps);
            P = zeros(Ntrials, numOfTimeSteps, model.getNumberOfSpecies());
        end
        T(n, :) = TArray;
        P(n, :, :) = XArray;
    end
    T = T(1, :);

    inputFunction = ch.ethz.bhepp.utils.SinusoidalFunction(input_offset, input_amplitude, input_frequency, initial_phase);
    model = ch.ethz.bhepp.ssasolver.models.LeloupGoldbeterCircadianClock(omega, inputFunction);
end
