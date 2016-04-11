function [T, X] = RunODE(t0, tf, recordStep, ...
    input_offset, input_amplitude, input_period, ...
    initial_phase)

    x0 = zeros(16, 1);

    % initial maximum
%     x0 = [4.5579, 2.8471, 7.5723, 1.1340, 7.6655, 0.1216, 0.6975, 4.6534, 1.6187, 0.2459, 0.2171, 2.0359, 0.9254, 0.8791, 0.3163, 0.6039];
    
%     % initial minimum
%     x0 = [0.0045, 0.0079, 8.5882, 0.2578, 3.0212, 0.0857, 0.7127, 2.1006, 3.4558, 0.2424, 0.2471, 1.6297, 0.8987, 0.2247, 0.1664, 0.2902];

    p = ModelParameters();

    input_frequency = 1 ./ input_period;
    inputFunction = @(t, x) input_offset ... 
        + input_amplitude * sin(2 * pi * t * input_frequency + initial_phase);

    [T, X] = ode15s(@ModelODE, t0:recordStep:tf, x0, odeset(), p, inputFunction);

    X = X(:, 1);

end
