function dx = CircadianClock(t, x, input_function)
    input = input_function(t, x);

%     um = 1e-6;
    um = 1.0;

    v_s = 0.76 * um;
    v_m = 0.65 * um;
    K_m = 0.5 * um;
    k_s = 0.38;
    v_d = 0.95 * um;
    k_1 = 1.9;
    k_2 = 1.3;
    K_I = 1 * um;
    K_d = 0.2 * um;
    n = 4;
    K_1 = 2 * um;
    K_2 = K_1;
    K_3 = K_1;
    K_4 = K_1;
    V_1 = 3.2 * um;
    V_2 = 1.58 * um;
    V_3 = 5 * um;
    V_4 = 2.5 * um;

    M = x(1);
    P_0 = x(2);
    P_1 = x(3);
    P_2 = x(4);
    P_N = x(5);

    term1 = V_1 * P_0 / (K_1 + P_0);
    term2 = V_2 * P_1 / (K_2 + P_1);
    term3 = V_3 * P_1 / (K_3 + P_1);
    term4 = V_4 * P_2 / (K_4 + P_2);
    term5 = k_1 * P_2;
    term6 = k_2 * P_N;

    dM = input * v_s * K_I^n / (K_I^n + P_N^n) - v_m * M / (K_m + M);
    dP_0 = k_s * M - term1 + term2;
    dP_1 = + term1 - term2 - term3 + term4;
    dP_2 = + term3 - term4 - term5 + term6 - v_d * P_2 / (K_d + P_2);
    dP_N = term5 - term6;

    dx = zeros(5, 1);
    dx(1) = dM;
    dx(2) = dP_0;
    dx(3) = dP_1;
    dx(4) = dP_2;
    dx(5) = dP_N;

end
