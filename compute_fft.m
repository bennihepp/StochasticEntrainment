function [omega, y] = compute_fft(x, dt)
    L = length(x);
    NFFT = pow2(nextpow2(L));
    y = fft(x, NFFT);

    Fs = 1 / dt;
    f = (0:NFFT-1) * (Fs / NFFT);
    f = f(1:NFFT/2+1);

    y = y(1:NFFT/2+1);
    omega = 2 * pi * f;
end
