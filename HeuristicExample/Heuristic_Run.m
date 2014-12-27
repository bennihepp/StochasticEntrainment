function [T, output] = Heuristic_Run(Ntrials, t0, tf, dt, A_0, A_1, epsilon, omega_0, omega_1, phi_0_rand, phi_1_rand)

    T = t0:dt:tf;
    output = zeros(length(T), Ntrials);
    for n=1:Ntrials
        phi_0 = phi_0_rand(n, 1, length(T));
        phi_1 = phi_1_rand(n, 1, length(T));
        output(:, n) = A_0 * cos(omega_0 * T + phi_0) + epsilon * A_1 * cos(omega_1 * T + phi_1);
    end

end
