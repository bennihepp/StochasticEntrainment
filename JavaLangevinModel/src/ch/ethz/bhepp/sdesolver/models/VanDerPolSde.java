package ch.ethz.bhepp.sdesolver.models;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class VanDerPolSde implements Sde {

	public final int NUM_OF_SPECIES = 2;
	public final int NUM_OF_REACTIONS = 4;

//	public final double Na = 6.022 * 1e23;

	private double B = 10.0;
	private double d = 2.0;

	private DoubleMatrix2D stoichiometryMatrix;
	private DoubleMatrix1D scaledPropensities;
	private double omega;
	private double sqrtOmega;
	private ScalarTrajectoryFunction inputFunction;


	public VanDerPolSde(double volume, double inputOffset, double inputAmplitude, double inputFrequency) {
		this(volume, new SinusoidalFunction(inputOffset, inputAmplitude, inputFrequency));
	}

	public VanDerPolSde(double volume, ScalarTrajectoryFunction inputFunction) {
//		omega = 1e-6 * Na * volume;
		omega = volume;
	    this.sqrtOmega = Math.sqrt(omega);
		this.inputFunction = inputFunction;
		stoichiometryMatrix = new DenseDoubleMatrix2D(NUM_OF_REACTIONS, NUM_OF_SPECIES);
		stoichiometryMatrix.set(0, 0, +1);
		stoichiometryMatrix.set(1, 1, -1);
		stoichiometryMatrix.set(2, 1, +1);
		stoichiometryMatrix.set(3, 1, -1);
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

		double x1 = X.get(0);
		double x2 = X.get(1);
		double input = inputFunction.evaluate(t, X);

		F.set(0, x2);
		F.set(1, input + (d - B * Math.pow(x1, 2.0)) * x2 - x1);

        scaledPropensities.set(0, Math.abs(x2));
        scaledPropensities.set(1, Math.abs(B * Math.pow(x1, 2.0) * x2));
        scaledPropensities.set(2, Math.abs(input + d * x2));
        scaledPropensities.set(3, Math.abs(x1));
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
