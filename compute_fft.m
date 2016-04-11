function [omega, y] = compute_fft(x, dt)
    L = length(x);
    np2 = nextpow2(L);
    if pow2(np2) > L
        np2 = np2 - 1;
    end
    NFFT = pow2(np2);
    y = fft(x, NFFT) / L;

    Fs = 1 / dt;
    f = Fs / 2 * linspace(0, 1, NFFT/2+1);

%     f = (0:NFFT-1) * (Fs / NFFT);
%     f = f(1:NFFT/2+1);

    y = y(1:NFFT/2+1);
    omega = 2 * pi * f;
end
