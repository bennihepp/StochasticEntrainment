function value = MultiplicativeForcing(t, x, period, amplitude)
    omega = 2 * pi / period;
    M = [0, 0; amplitude, 0];
    value = (M * x) * sin(omega * t);
end
