function [T, output] = VanDerPol_Run_Java(Ntrials, t0, tf, dt, recordStep, ...
    volume, input_offset, input_amplitude, input_period, initial_phase)

    if nargin < 10
        initial_phase = 0;
    end

    x0 = [1.0; 1.0];

    input_frequency = 1 ./ input_period;

    do_parallel = true;
%     do_parallel = false;
%     if Ntrials > 1
%         do_parallel = true;
%     end

    if do_parallel
        [T, X, ~] = SolveS_Java_Parallel(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, initial_phase, Ntrials, ...
            recordStep);
    else
        [T, X, ~] = SolveS_Java(x0, tf, dt, volume, ...
            input_offset, input_amplitude, input_frequency, initial_phase, Ntrials, ...
            recordStep);
    end


    output = squeeze(X(:, :, 1));
    output = transpose(output);

end
