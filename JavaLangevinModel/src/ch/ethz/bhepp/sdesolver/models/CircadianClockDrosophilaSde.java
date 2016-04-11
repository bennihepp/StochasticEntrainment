package ch.ethz.bhepp.sdesolver.models;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class CircadianClockDrosophilaSde implements Sde {

	public final int NUM_OF_SPECIES = 5;
	public final int NUM_OF_REACTIONS = 10;

	public final double Na = 6.022 * 1e23;

//	private double um = 1.0e-6;
	private double um = 1.0;
	private double v_s = 0.76 * um;
	private double v_m = 0.65 * um;
	private double K_m = 0.5 * um;
	private double k_s = 0.38;
	private double v_d = 0.95 * um;
	private double k_1 = 1.9;
	private double k_2 = 1.3;
	private double K_I = 1.0 * um;
	private double K_d = 0.2 * um;
	private double n = 4;
	private double K_1 = 2.0 * um;
	private double K_2 = K_1;
	private double K_3 = K_1;
	private double K_4 = K_1;
	private double V_1 = 3.2 * um;
	private double V_2 = 1.58 * um;
	private double V_3 = 5.0 * um;
	private double V_4 = 2.5 * um;

	private DoubleMatrix2D stoichiometryMatrix;
	private DoubleMatrix1D scaledPropensities;
	private double omega;
	private double sqrtOmega;
	private ScalarTrajectoryFunction inputFunction;


	public CircadianClockDrosophilaSde(double volume, double inputOffset, double inputAmplitude, double inputFrequency) {
		this(volume, new SinusoidalFunction(inputOffset, inputAmplitude, inputFrequency));
	}

	public CircadianClockDrosophilaSde(double volume, ScalarTrajectoryFunction inputFunction) {
		omega = um * Na * volume;
	    this.sqrtOmega = Math.sqrt(omega);
		this.inputFunction = inputFunction;
		stoichiometryMatrix = new DenseDoubleMatrix2D(NUM_OF_REACTIONS, NUM_OF_SPECIES);
		stoichiometryMatrix.set(0, 0, +1); // term1
		stoichiometryMatrix.set(1, 0, -1); // term2
		stoichiometryMatrix.set(2, 1, +1); // term3
		stoichiometryMatrix.set(3, 1, -1); // term4
		stoichiometryMatrix.set(3, 2, +1); // term4
		stoichiometryMatrix.set(4, 1, +1); // term5
		stoichiometryMatrix.set(4, 2, -1); // term5
		stoichiometryMatrix.set(5, 2, -1); // term6
		stoichiometryMatrix.set(5, 3, +1); // term6
		stoichiometryMatrix.set(6, 2, +1); // term7
		stoichiometryMatrix.set(6, 3, -1); // term7
		stoichiometryMatrix.set(7, 3, -1); // term8
		stoichiometryMatrix.set(7, 4, +1); // term8
		stoichiometryMatrix.set(8, 3, +1); // term9
		stoichiometryMatrix.set(8, 4, -1); // term9
		stoichiometryMatrix.set(9, 3, -1); // term10
		scaledPropensities = new DenseDoubleMatrix1D(NUM_OF_REACTIONS);
	}

	public int getDriftDimension() {
		return NUM_OF_SPECIES;
	}

	public int getDiffusionDimension() {
		return NUM_OF_REACTIONS;
	}

	public double getOmega() {
		return omega;
	}

	public void computeDriftAndDiffusion(double t, DoubleMatrix1D X,
			DoubleMatrix1D F, DoubleMatrix2D G) throws Exception {

		double M = X.get(0);
		double P_0 = X.get(1);
		double P_1 = X.get(2);
		double P_2 = X.get(3);
		double P_N = X.get(4);
		double input = inputFunction.evaluate(t, X);

		double term1 = input * v_s * Math.pow(K_I, n) / (Math.pow(K_I, n) + Math.pow(P_N, n));
		double term2 = v_m * M / (K_m + M);
		double term3 = k_s * M;
		double term4 = V_1 * P_0 / (K_1 + P_0);
		double term5 = V_2 * P_1 / (K_2 + P_1);
		double term6 = V_3 * P_1 / (K_3 + P_1);
		double term7 = V_4 * P_2 / (K_4 + P_2);
		double term8 = k_1 * P_2;
		double term9 = k_2 * P_N;
		double term10 = v_d * P_2 / (K_d + P_2);

		scaledPropensities.set(0, term1);
		scaledPropensities.set(1, term2);
		scaledPropensities.set(2, term3);
		scaledPropensities.set(3, term4);
		scaledPropensities.set(4, term5);
		scaledPropensities.set(5, term6);
		scaledPropensities.set(6, term7);
		scaledPropensities.set(7, term8);
		scaledPropensities.set(8, term9);
		scaledPropensities.set(9, term10);

		F.assign(0.0);
		G.assign(0.0);
		for (int k=0; k < NUM_OF_REACTIONS; k++) {
			double propensity = scaledPropensities.get(k);
			double propensitySqrt = Math.sqrt(propensity);
			for (int i=0; i < NUM_OF_SPECIES; i++) {
				double newValueF = F.get(i) + propensity * stoichiometryMatrix.get(k, i);
				F.set(i, newValueF);
				double newValue = G.get(i, k) + propensitySqrt * stoichiometryMatrix.get(k, i) / sqrtOmega;
				G.set(i, k, newValue);
			}
		}

//		F.set(0, term1);
//		F.set(1, term2);
//		F.set(2, term3);
//		F.set(3, term4);
//		F.set(4, term5);
//		F.set(5, term6);
//		F.set(6, term7);
//		F.set(7, term8);
//		F.set(8, term9);
//		F.set(9, term10);

//		G.assign(0.0);
//		for (int k=0; k < NUM_OF_REACTIONS; k++) {
//			double propensitySqrt = Math.sqrt(scaledPropensities.get(k));
//			for (int i=0; i < NUM_OF_SPECIES; i++) {
//				double newValue = G.get(i, k) + propensitySqrt * stoichiometryMatrix.get(k, i) / sqrtOmega;
//				G.set(i, k, newValue);
//			}
//		}

	}

}
