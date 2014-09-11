function [T, output] = VanDerPol_Run2(Ntrials, t0, tf, dt, recordStep, volume, input_offset, input_amplitude, input_period)

    x0 = [1.0; 1.0];

    input_frequency = 1 ./ input_period;

    do_parallel = false;
    if Ntrials > 1
        do_parallel = true;
    end

    if do_parallel
        [T, X, ~] = SolveS_Java_Parallel2(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, Ntrials, ...
            recordStep);
    else
        [T, X, ~] = SolveS_Java(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, Ntrials, ...
            recordStep);
    end


    output = squeeze(X(:, :, 1));
    output = transpose(output);

end
