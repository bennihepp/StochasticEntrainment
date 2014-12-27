package ch.ethz.bhepp.sdesolver.models;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class NFkBSpikySde implements Sde {

	public final int NUM_OF_SPECIES = 3;
	public final int NUM_OF_REACTIONS = 5;

	public final double Na = 6.022 * 1e23;

	private double A = 0.007;
	private double B = 954.5;

	private double C0 = 0.035;
	private double delta = 0.029;
	private double epsilon = 2 * 1e-5;

	private DoubleMatrix2D stoichiometryMatrix;
	private DoubleMatrix1D scaledPropensities;
	private double omega;
	private double sqrtOmega;
	private ScalarTrajectoryFunction inputFunction;


	public NFkBSpikySde(double volume, double inputOffset, double inputAmplitude, double inputFrequency) {
		this(volume, new SinusoidalFunction(inputOffset, inputAmplitude, inputFrequency));
	}

	public NFkBSpikySde(double volume, ScalarTrajectoryFunction inputFunction) {
		omega = 1e-6 * Na * volume;
	    this.sqrtOmega = Math.sqrt(omega);
		this.inputFunction = inputFunction;
		stoichiometryMatrix = new DenseDoubleMatrix2D(NUM_OF_REACTIONS, NUM_OF_SPECIES);
		stoichiometryMatrix.set(0, 0, +1);
		stoichiometryMatrix.set(1, 0, -1);
		stoichiometryMatrix.set(2, 1, +1);
		stoichiometryMatrix.set(3, 1, -1);
		stoichiometryMatrix.set(3, 2, +1);
		stoichiometryMatrix.set(4, 2, -1);
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

		double Nn = X.get(0);
		double Im = X.get(1);
		double I = X.get(2);
		double C = C0 * inputFunction.evaluate(t, X);

		F.set(0, A * (1 - Nn) / (epsilon + I) - B * (I * Nn) / (delta + Nn));
		F.set(1, Nn * Nn - Im);
		F.set(2, Im - C * (1 - Nn) * I / (epsilon + I));

        scaledPropensities.set(0, A * (1 - Nn) / (epsilon + I));
        scaledPropensities.set(1, B * (I * Nn) / (delta + Nn));
        scaledPropensities.set(2, Nn * Nn);
        scaledPropensities.set(3, Im);
        scaledPropensities.set(4, C * (1 - Nn) * I / (epsilon + I));
//        for (int k=0; k < scaledPropensities.size(); k++)
//        	if (scaledPropensities.get(k) < 0.0)
//        		scaledPropensities.set(k, 0.0);

		G.assign(0.0);
		for (int k=0; k < NUM_OF_REACTIONS; k++) {
			double propensitySqrt = Math.sqrt(scaledPropensities.get(k));
			for (int i=0; i < NUM_OF_SPECIES; i++) {
				double newValue = G.get(i, k) + propensitySqrt * stoichiometryMatrix.get(k, i) / sqrtOmega;
				G.set(i, k, newValue);
			}
		}

	}

}
