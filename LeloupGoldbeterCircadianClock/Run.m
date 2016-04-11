function [T, output, model] = Run(Ntrials, t0, tf, recordStep, omega, ...
    input_offset, input_amplitude, input_period, initial_phase, ...
    printMessages)

    if nargin < 10
        printMessages = false;
    end

%     x0 = [2.7438, 1.7099, 4.5363, 0.6806, 4.5912, 0.0731, 0.4154, 2.7857, 0.9644, 0.1472, 0.1306, 1.2233, 0.5518, 0.5328, 0.1914, 0.3622]';
    x0 = zeros(16, 1);

    input_frequency = 1 ./ input_period;

%     do_parallel = false;
%     if Ntrials > 1
%         do_parallel = true;
%     end
    do_parallel = 1;

    if do_parallel
        [T, X, model] = SolveS_Java_Parallel(x0, t0, tf, recordStep, ...
                            omega, input_offset, input_amplitude, ...
                            input_frequency, initial_phase, Ntrials, ...
                            printMessages);
    else
        [T, X, model] = SolveS_Java(x0, t0, tf, recordStep, omega, ...
            input_offset, input_amplitude, input_frequency, initial_phase, ...
            Ntrials);
    end

    output = X';

%     per_indices = [0] + 1;
%     output = sum(X(:, :, per_indices), 3)';
%     output = squeeze(sum(X(:, :, per_indices), 3));
%     output = squeeze(X);

%     output = squeeze(sum(X(:, :, 2:5), 3));
%     output = transpose(output);

end
