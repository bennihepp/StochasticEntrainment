function value = AdditiveForcing(t, x, period, amplitude)
    omega = 2 * pi / period;
    value = [0; amplitude] * sin(omega * t);
end
