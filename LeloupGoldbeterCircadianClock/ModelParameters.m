function [p] = Parameters()

    p = struct();

	p.k_1 = 0.4;
	p.k_2 = 0.2;
	p.k_3 = 0.4;
	p.k_4 = 0.2;
	p.k_5 = 0.4;
	p.k_6 = 0.2;
	p.k_7 = 0.5;
	p.k_8 = 0.1;
	p.K_AP = 0.7;
	p.K_AC = 0.6;
	p.K_IB = 2.2;
	p.k_dmb = 0.01;
	p.k_dmc = 0.01;
	p.k_dmp = 0.01;
	p.k_dn = 0.01;
	p.k_dnc = 0.12;
	p.K_d = 0.3;
	p.K_dp = 0.1;
	p.K_p = 0.1;
	p.K_mB = 0.4;
	p.K_mC = 0.4;
	p.K_mP = 0.31;
	p.k_sB = 0.12;
	p.k_sC = 1.6;
	p.k_sP = 0.6;
	p.m = 2;
	p.n = 4;
	p.V_1B = 0.5;
	p.V_1C = 0.6;
	p.V_1P = 0.4;
	p.V_1PC = 0.4;
	p.V_2B = 0.1;
	p.V_2C = 0.1;
	p.V_2P = 0.3;
	p.V_2PC = 0.1;
	p.V_3B = 0.5;
	p.V_3PC = 0.4;
	p.V_4B = 0.2;
	p.V_4PC = 0.1;
	p.V_phos = 0.4;
	p.v_dBC = 0.5;
	p.v_dBN = 0.6;
	p.v_dCC = 0.7;
	p.v_dIN = 0.8;
	p.v_dPC = 0.7;
	p.v_dPCC = 0.7;
	p.v_dPCN = 0.7;
	p.v_mB = 0.8;
	p.v_mC = 1.0;
	p.v_mP = 1.1;
	p.v_sB = 1.0;
	p.v_sC = 1.1;
	p.v_sP = 1.5;

	p.M_P = 0 + 1;
	p.M_C = 1 + 1;
	p.M_B = 2 + 1;
	p.P_C = 3 + 1;
	p.C_C = 4 + 1;
	p.P_CP = 5 + 1;
	p.C_CP = 6 + 1;
	p.PC_C = 7 + 1;
	p.PC_N = 8 + 1;
	p.PC_CP = 9 + 1;
	p.PC_NP = 10 + 1;
	p.B_C = 11 + 1;
	p.B_CP = 12 + 1;
	p.B_N = 13 + 1;
	p.B_NP = 14 + 1;
	p.I_N = 15 + 1;

end