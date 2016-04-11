function [T, output, omega, sde] = Run(Ntrials, t0, tf, dt, recordStep, ...
    volume, input_offset, input_amplitude, input_period, ...
    initial_phase, printMessages)

    if nargin < 10
        initial_phase = 0;
    end

    if nargin < 11
        printMessages = false;
    end

%     x0 = [4.5579, 2.8471, 7.5723, 1.1340, 7.6655, 0.1216, 0.6975, 4.6534, 1.6187, 0.2459, 0.2171, 2.0359, 0.9254, 0.8791, 0.3163, 0.6039]';
    x0 = 1e3 * [2.7438, 1.7099, 4.5363, 0.6806, 4.5912, 0.0731, 0.4154, 2.7857, 0.9644, 0.1472, 0.1306, 1.2233, 0.5518, 0.5328, 0.1914, 0.3622]';

    input_frequency = 1 ./ input_period;

    do_parallel = false;
    if Ntrials > 1
        do_parallel = true;
    end

    if do_parallel
        [T, X] = SolveS_Java_Parallel(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, initial_phase, ...
            Ntrials, recordStep, printMessages);
    else
        [T, X, omega, sde] = SolveS_Java(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, initial_phase, ...
            Ntrials, recordStep);
    end

%     per_indices = [8, 10] + 1;
%     output = squeeze(sum(X(:, :, per_indices), 3));
    output = squeeze(X);

%     output = squeeze(sum(X(:, :, 2:5), 3));
%     output = transpose(output);

end
