function [T, output] = VanDerPol_Run(Ntrials, t0, tf, dt, omega, ...
    additive_forcing_func, multiplicative_forcing_func)

    y0 = [1; 1];

    SDE = VanDerPol_Model(y0, omega, additive_forcing_func, multiplicative_forcing_func);
    Nsteps = ceil((tf - t0) / dt);

    [Paths, Times, Z] = SDE.simByEuler(Nsteps, 'DeltaTime', dt, 'NTRIALS', Ntrials);

    T = Times;
    output = squeeze(Paths(:, 1, :));

end
