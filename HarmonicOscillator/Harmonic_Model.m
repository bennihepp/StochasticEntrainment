function SDE = Harmonic_Model(x0, omega, additive_input_func, multiplicative_input_func)

    [k, c] = Parameters();

    function drift = F(t, y)

        additive_input = additive_input_func(t, y);
        multiplicative_input = multiplicative_input_func(t, y);

        dy = zeros(2, 1);

        dy(1) = y(2);
        dy(2) = - k * y(1) - c * y(2);

        drift = dy + additive_input + multiplicative_input;

    end

    function diff = G(~, ~)
        diff = eye(2);
        diff = diff ./ sqrt(omega);
    end

    SDE = sde(@F, @G, 'StartState', x0);

end
