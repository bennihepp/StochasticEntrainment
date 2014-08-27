function [T, output] = CircadianClock_Run(Ntrials, t0, tf, dt, recordStep, volume, input_offset, input_amplitude, input_period)

    um = 1.0;
    x0 = [0.1; 0.25; 0.25; 0.25; 0.25] * um;

    input_frequency = 1 ./ input_period;

    if Ntrials > 1
        [T, X, ~] = SolveS_Java_Parallel(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, Ntrials, ...
            recordStep);
    else
        [T, X, ~] = SolveS_Java(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, Ntrials, ...
            recordStep);
    end


    output = squeeze(sum(X(:, :, 2:5), 3));
    output = transpose(output);

end
