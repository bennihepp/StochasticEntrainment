function [omega, y] = compute_normalized_fft_truncated(x, dt, min_omega, max_omega)
    [omega, y] = compute_normalized_fft(x, dt);
    i1 = find(omega < min_omega, 1, 'last');
    i2 = find(omega > max_omega, 1, 'first');
    if isempty(i1)
        i1 = 1;
    end
    if isempty(i2)
        i2 = length(omega);
    end
    omega = omega(i1:i2);
    y = y(i1:i2);
end
