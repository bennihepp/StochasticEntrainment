% volume in liters
function SDE = SpikyModelS(x0, omega, input_function)

    %% General definitions
    A = 0.007;
    B = 954.5;
%     C = 0.035;
    C0 = 0.035;
    delta = 0.029;
    epsilon = 2 * 1e-5;

    function prop = scaled_prop_func(t, x)
        input = input_function(t, x);
        C = C0 * input;
        prop = zeros(5, 1);
        prop(1) = A * (1 - Nn) / (epsilon + I);
        prop(2) = B * (I * Nn) / (delta + Nn);
        prop(3) = Nn.^2;
        prop(4) = Im;
        prop(5) = C * (1 - Nn) * I / (epsilon + I);
    end

    stoch_matrix = zeros(5, 3);
    stoch_matrix(1, 1) = +1;
    stoch_matrix(2, 1) = -1;
    stoch_matrix(3, 2) = +1;
    stoch_matrix(4, 2) = -1;
    stoch_matrix(4, 3) = +1;
    stoch_matrix(5, 3) = -1;

    F = @(t, x) SpikyModel(t, x, input_function);

    function diff = G(t, x)
        diff = zeros([size(stoch_matrix, 2), size(stoch_matrix, 1)]);
        propensities = scaled_prop_func(t, x);
        propensities(propensities < 0) = 0;
        for k=1:size(diff, 2)
            diff(:, k) = diff(:, k) + stoch_matrix(k, :)' * sqrt(propensities(k));
        end
        diff = diff ./ sqrt(omega);
    end

    SDE = sde(F, @G, 'StartState', x0);

end
