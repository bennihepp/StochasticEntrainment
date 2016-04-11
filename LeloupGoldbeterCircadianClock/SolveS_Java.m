% volume in liters
% time in minutes
function [T, P, model] = SolveS_Java(x0, t0, tf, recordStep, omega, ...
                            input_offset, input_amplitude, ...
                            input_frequency, initial_phase, Ntrials, ...
                            seed)

    if nargin < 10
        Ntrials = 1;
    end

    if nargin < 11
        seed = randi([0, 2.^31-2]);
%         seed = 4942634935159178;
%         disp(['Java seed=', int2str(seed)]);
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
    model = ch.ethz.bhepp.ssasolver.models.LeloupGoldbeterCircadianClock(omega, inputFunction);

    stepper = ch.ethz.bhepp.ssasolver.SimpleSSAStepper(model, seed);
    solver = ch.ethz.bhepp.ssasolver.SSASolver(stepper);

    x0Matrix = ch.ethz.bhepp.utils.MatrixHelper.createZeroDoubleMatrix1D(size(x0, 1));
    if size(x0, 2) == 1
        x0 = repmat(x0, 1, Ntrials);
    end

    for i=1:size(x0, 1)
        x0Matrix.set(i-1, x0(i, 1));
    end
    solution = solver.solve(t0, x0Matrix, tf, recordStep);
    TArray = solution.getTArray();
    XArray = solution.getXArray();
%     T = TArray;
%     P = zeros(Ntrials, size(XArray, 1), size(XArray, 2));
    numOfTimeSteps = length(TArray);
    T = zeros(Ntrials, numOfTimeSteps);
    P = zeros(Ntrials, numOfTimeSteps, model.getNumberOfSpecies());

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
