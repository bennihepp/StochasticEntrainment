package ch.ethz.bhepp.ssasolver.models;

import ch.ethz.bhepp.ssasolver.StochasticReactionNetworkModel;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class LeloupGoldbeterCircadianClock implements StochasticReactionNetworkModel {

	public final int NUM_OF_SPECIES = 16;
	public final int NUM_OF_REACTIONS = 52;

	private final static LeloupGoldbeterCircadianClockParameters p
		= new LeloupGoldbeterCircadianClockParameters();

	private ScalarTrajectoryFunction inputFunction;
	private double omega;

	public LeloupGoldbeterCircadianClock(double omega, double inputOffset,
			double inputAmplitude, double inputFrequency) {
		this(omega, new SinusoidalFunction(inputOffset, inputAmplitude, inputFrequency));
	}

	public LeloupGoldbeterCircadianClock(double omega, ScalarTrajectoryFunction inputFunction) {
		this.omega = omega;
		this.inputFunction = inputFunction;
	}

	public LeloupGoldbeterCircadianClockParameters getParameters() {
		return p;
	}

	@Override
	public int getNumberOfSpecies() {
		return NUM_OF_SPECIES;
	}

	@Override
	public int getNumberOfReactions() {
		return NUM_OF_REACTIONS;
	}

	@Override
	public double computePropensity(int reaction, double t, double[] x) {
		if (reaction == 0) {
			double B_N_pow = Math.pow(x[p.B_N], p.n);
			double v = (p.v_sP * omega) * B_N_pow / (Math.pow(p.K_AP * omega, p.n) + B_N_pow);
			return v * inputFunction.evaluate(t, x);
		} else if (reaction == 1) {
			double M_P = x[p.M_P];
			return (p.v_mP * omega) * M_P / (p.K_mP * omega + M_P);
		} else if (reaction == 2) {
			double M_P = x[p.M_P];
			return p.k_dmp * M_P;
		} else if (reaction == 3) {
			double B_N_pow = Math.pow(x[p.B_N], p.n);
			return (p.v_sC * omega)* B_N_pow / (Math.pow(p.K_AC * omega, p.n) + B_N_pow);
		} else if (reaction == 4) {
			double M_C = x[p.M_C];
			return (p.v_mC * omega) * M_C / (p.K_mC * omega + M_C);
		} else if (reaction == 5) {
			double M_C = x[p.M_C];
			return p.k_dmc * M_C;
		} else if (reaction == 6) {
			double B_N_pow = Math.pow(x[p.B_N], p.m);
			double K_IB_pow = Math.pow(p.K_IB * omega, p.m);
			return (p.v_sB * omega) * K_IB_pow / (K_IB_pow + B_N_pow);
		} else if (reaction == 7) {
			double M_B = x[p.M_B];
			return (p.v_mB * omega)* M_B / (p.K_mB * omega + M_B); 
		} else if (reaction == 8) {
			double M_B = x[p.M_B];
			return p.k_dmb * M_B;
		} else if (reaction == 9) {
			double M_P = x[p.M_P];
			return p.k_sP * M_P;
		} else if (reaction == 10) {
			double P_C = x[p.P_C];
			return (p.V_1P * omega) * P_C / (p.K_p * omega + P_C);
		} else if (reaction == 11) {
			double P_CP = x[p.P_CP];
			return (p.V_2P * omega) * P_CP / (p.K_dp * omega + P_CP);
		} else if (reaction == 12) {
			double PC_C = x[p.PC_C];
			return p.k_4 * PC_C;
		} else if (reaction == 13) {
			double P_C = x[p.P_C];
			double C_C = x[p.C_C];
			return (p.k_3 / omega) * P_C * C_C;
		} else if (reaction == 14) {
			double P_C = x[p.P_C];
			return p.k_dn * P_C;
		} else if (reaction == 15) {
			double M_C = x[p.M_C];
			return p.k_sC * M_C;
		} else if (reaction == 16) {
			double C_C = x[p.C_C];
			return (p.V_1C * omega) * C_C / (p.K_p * omega + C_C);
		} else if (reaction == 17) {
			double C_CP = x[p.C_CP];
			return (p.V_2C * omega) * C_CP / (p.K_dp * omega + C_CP);
		} else if (reaction == 18) {
			double C_C = x[p.C_C];
			return p.k_dnc * C_C;
		} else if (reaction == 19) {
			double P_CP = x[p.P_CP];
			return (p.v_dPC * omega) * P_CP / (p.K_d * omega + P_CP);
		} else if (reaction == 20) {
			double P_CP = x[p.P_CP];
			return p.k_dn * P_CP;
		} else if (reaction == 21) {
			double C_CP = x[p.C_CP];
			return (p.v_dCC * omega) * C_CP / (p.K_d * omega + C_CP);
		} else if (reaction == 22) {
			double C_CP = x[p.C_CP];
			return p.k_dn * C_CP;
		} else if (reaction == 23) {
			double PC_C = x[p.PC_C];
			return (p.V_1PC * omega) * PC_C / (p.K_p * omega + PC_C);
		} else if (reaction == 24) {
			double PC_CP = x[p.PC_CP];
			return (p.V_2PC * omega) * PC_CP / (p.K_dp * omega + PC_CP);
		} else if (reaction == 25) {
			double PC_N = x[p.PC_N];
			return p.k_2 * PC_N;
		} else if (reaction == 26) {
			double PC_C = x[p.PC_C];
			return p.k_1 * PC_C;
		} else if (reaction == 27) {
			double PC_C = x[p.PC_C];
			return p.k_dn * PC_C;
		} else if (reaction == 28) {
			double PC_N = x[p.PC_N];
			return (p.V_3PC * omega) * PC_N / (p.K_p * omega + PC_N);
		} else if (reaction == 29) {
			double PC_NP = x[p.PC_NP];
			return (p.V_4PC * omega) * PC_NP / (p.K_dp * omega + PC_NP);
		} else if (reaction == 30) {
			double PC_N = x[p.PC_N];
			double B_N = x[p.B_N];
			return (p.k_7 / omega)* PC_N * B_N; 
		} else if (reaction == 31) {
			double I_N = x[p.I_N];
			return p.k_8 * I_N;
		} else if (reaction == 32) {
			double PC_N = x[p.PC_N];
			return p.k_dn * PC_N;
		} else if (reaction == 33) {
			double PC_CP = x[p.PC_CP];
			return (p.V_dPCC * omega) * PC_CP / (p.K_d * omega + PC_CP);
		} else if (reaction == 34) {
			double PC_CP = x[p.PC_CP];
			return p.k_dn * PC_CP;
		} else if (reaction == 35) {
			double PC_NP = x[p.PC_NP];
			return (p.V_dPCN * omega) * PC_NP / (p.K_d * omega + PC_NP);
		} else if (reaction == 36) {
			double PC_NP = x[p.PC_NP];
			return p.k_dn * PC_NP;
		} else if (reaction == 37) {
			double M_B = x[p.M_B];
			return p.k_sB * M_B;
		} else if (reaction == 38) {
			double B_C = x[p.B_C];
			return (p.V_1B * omega) * B_C / (p.K_p * omega + B_C);
		} else if (reaction == 39) {
			double B_CP = x[p.B_CP];
			return (p.V_2B * omega) * B_CP / (p.K_dp * omega + B_CP);
		} else if (reaction == 40) {
			double B_C = x[p.B_C];
			return p.k_5 * B_C;
		} else if (reaction == 41) {
			double B_N = x[p.B_N];
			return p.k_6 * B_N;
		} else if (reaction == 42) {
			double B_C = x[p.B_C];
			return p.k_dn * B_C;
		} else if (reaction == 43) {
			double B_CP = x[p.B_CP];
			return (p.V_dBC * omega) * B_CP / (p.K_d * omega + B_CP);
		} else if (reaction == 44) {
			double B_CP = x[p.B_CP];
			return p.k_dn * B_CP;
		} else if (reaction == 45) {
			double B_N = x[p.B_N];
			return (p.V_3B * omega) * B_N / (p.K_p * omega + B_N);
		} else if (reaction == 46) {
			double B_NP = x[p.B_NP];
			return (p.V_4B * omega) * B_NP / (p.K_dp * omega + B_NP);
		} else if (reaction == 47) {
			double B_N = x[p.B_N];
			return p.k_dn * B_N;
		} else if (reaction == 48) {
			double B_NP = x[p.B_NP];
			return (p.V_dBN * omega) * B_NP / (p.K_d * omega + B_NP);
		} else if (reaction == 49) {
			double B_NP = x[p.B_NP];
			return p.k_dn * B_NP;
		} else if (reaction == 50) {
			double I_N = x[p.I_N];
			return (p.V_dIN * omega) * I_N / (p.K_d * omega + I_N);
		} else if (reaction == 51) {
			double I_N = x[p.I_N];
			return p.k_dn * I_N;
//		} else if (reaction == 52) {
//		} else if (reaction == 53) {
		} else {
			throw new IllegalArgumentException();
		}
	}

	public double[] computePropensities(double t, double[] x) {
		double[] propensities = new double[getNumberOfReactions()];
		computePropensities(t, x, propensities);
		return propensities;
	}

	@Override
	public void computePropensities(double t, double[] x, double[] propensities) {
		computePropensitiesAndSum(t, x, propensities);
	}

	@Override
	public double computePropensitiesAndSum(double t, double[] x, double[] propensities) {
		double propSum = 0.0;
		for (int r=0; r < getNumberOfReactions(); r++) {
			propensities[r] = computePropensity(r, t, x);
			propSum += propensities[r];
		}
		return propSum;
	}

	@Override
	public void changeState(int reaction, double t, double[] x) {
		if (reaction == 0) {
			x[p.M_P] += 1;
		} else if (reaction == 1) {
			x[p.M_P] -= 1;
		} else if (reaction == 2) {
			x[p.M_P] -= 1;
		} else if (reaction == 3) {
			x[p.M_C] += 1;
		} else if (reaction == 4) {
			x[p.M_C] -= 1;
		} else if (reaction == 5) {
			x[p.M_C] -= 1;
		} else if (reaction == 6) {
			x[p.M_B] += 1;
		} else if (reaction == 7) {
			x[p.M_B] -= 1;
		} else if (reaction == 8) {
			x[p.M_B] -= 1;
		} else if (reaction == 9) {
			x[p.P_C] += 1;
		} else if (reaction == 10) {
			x[p.P_C] -= 1;
			x[p.P_CP] += 1;
		} else if (reaction == 11) {
			x[p.P_C] += 1;
			x[p.P_CP] -= 1;
		} else if (reaction == 12) {
			x[p.P_C] += 1;
			x[p.C_C] += 1;
			x[p.PC_C] -= 1;
		} else if (reaction == 13) {
			x[p.P_C] -= 1;
			x[p.C_C] -= 1;
			x[p.PC_C] += 1;
		} else if (reaction == 14) {
			x[p.P_C] -= 1;
		} else if (reaction == 15) {
			x[p.C_C] += 1;
		} else if (reaction == 16) {
			x[p.C_C] -= 1;
			x[p.C_CP] += 1;
		} else if (reaction == 17) {
			x[p.C_C] += 1;
			x[p.C_CP] -= 1;
		} else if (reaction == 18) {
			x[p.C_C] -= 1;
		} else if (reaction == 19) {
			x[p.P_CP] -= 1;
		} else if (reaction == 20) {
			x[p.P_CP] -= 1;
		} else if (reaction == 21) {
			x[p.C_CP] -= 1;
		} else if (reaction == 22) {
			x[p.C_CP] -= 1;
		} else if (reaction == 23) {
			x[p.PC_C] -= 1;
			x[p.PC_CP] += 1;
		} else if (reaction == 24) {
			x[p.PC_C] += 1;
			x[p.PC_CP] -= 1;
		} else if (reaction == 25) {
			x[p.PC_C] += 1;
			x[p.PC_N] -= 1;
		} else if (reaction == 26) {
			x[p.PC_C] -= 1;
			x[p.PC_N] += 1;
		} else if (reaction == 27) {
			x[p.PC_C] -= 1;
		} else if (reaction == 28) {
			x[p.PC_N] -= 1;
			x[p.PC_NP] += 1;
		} else if (reaction == 29) {
			x[p.PC_N] += 1;
			x[p.PC_NP] -= 1;
		} else if (reaction == 30) {
			x[p.PC_N] -= 1;
			x[p.B_N] -= 1;
			x[p.I_N] += 1;
		} else if (reaction == 31) {
			x[p.PC_N] += 1;
			x[p.B_N] += 1;
			x[p.I_N] -= 1;
		} else if (reaction == 32) {
			x[p.PC_N] -= 1;
		} else if (reaction == 33) {
			x[p.PC_CP] -= 1;
		} else if (reaction == 34) {
			x[p.PC_CP] -= 1;
		} else if (reaction == 35) {
			x[p.PC_NP] -= 1;
		} else if (reaction == 36) {
			x[p.PC_NP] -= 1;
		} else if (reaction == 37) {
			x[p.B_C] += 1;
		} else if (reaction == 38) {
			x[p.B_C] -= 1;
			x[p.B_CP] += 1;
		} else if (reaction == 39) {
			x[p.B_C] += 1;
			x[p.B_CP] -= 1;
		} else if (reaction == 40) {
			x[p.B_C] -= 1;
			x[p.B_N] += 1;
		} else if (reaction == 41) {
			x[p.B_C] += 1;
			x[p.B_N] -= 1;
		} else if (reaction == 42) {
			x[p.B_C] -= 1;
		} else if (reaction == 43) {
			x[p.B_CP] -= 1;
		} else if (reaction == 44) {
			x[p.B_CP] -= 1;
		} else if (reaction == 45) {
			x[p.B_N] -= 1;
			x[p.B_NP] += 1;
		} else if (reaction == 46) {
			x[p.B_N] += 1;
			x[p.B_NP] -= 1;
		} else if (reaction == 47) {
			x[p.B_N] -= 1;
		} else if (reaction == 48) {
			x[p.B_NP] -= 1;
		} else if (reaction == 49) {
			x[p.B_NP] -= 1;
		} else if (reaction == 50) {
			x[p.I_N] -= 1;
		} else if (reaction == 51) {
			x[p.I_N] -= 1;
//		} else if (reaction == 52) {
//			x[p.CB_star] += 1;
//		} else if (reaction == 53) {
//			x[p.CB_star] -= 1;
		} else {
			throw new IllegalArgumentException();
		}
	}

}
