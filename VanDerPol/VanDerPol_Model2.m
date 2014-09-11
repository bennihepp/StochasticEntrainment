function SDE = VanDerPol_Model2(x0, omega, additive_input_func, multiplicative_input_func)

    [B, d] = Parameters();

    function prop = scaled_prop_func(t, y, TNF)
        prop = zeros(4, 1);
        prop(1) = abs(y(2));
        prop(2) = abs((B * y(1).^2) * y(2));
        prop(3) = abs((d) * y(2));
        prop(4) = abs(y(1));
    end

    stoch_matrix = zeros(4, 2);
    stoch_matrix(1, 1) = +1;
    stoch_matrix(2, 2) = -1;
    stoch_matrix(3, 2) = +1;
    stoch_matrix(4, 2) = -1;

    function drift = F(t, y)

        additive_input = additive_input_func(t, y);
        multiplicative_input = multiplicative_input_func(t, y);

        dy = zeros(2, 1);

        propensities = scaled_prop_func(t, x);
        propensities(3) = propensities(3) + additive_input(2) + multiplicative_input(2);
        propensities(propensities < 0) = 0;

        for k=1:size(diff, 2)
            diff(:, k) = diff(:, k) + stoch_matrix(k, :)' * sqrt(propensities(k));
        end

%         dy(1) = y(2);
%         dy(2) = -(B * y(1).^2 - d) * y(2) - y(1);

%         drift = dy + additive_input + multiplicative_input;

    end

    function diff = G(t, x)
%         diff = sigma*eye(5);
        diff = zeros([size(stoch_matrix, 2), size(stoch_matrix, 1)]);
        propensities = scaled_prop_func(t, x);
        propensities(propensities < 0) = 0;
        for k=1:size(diff, 2)
            diff(:, k) = diff(:, k) + stoch_matrix(k, :)' * sqrt(propensities(k));
        end
        diff = diff ./ sqrt(omega);
    end

    SDE = sde(@F, @G, 'StartState', x0);

end
