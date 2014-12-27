package ch.ethz.bhepp.sdesolver.models;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.ssasolver.models.LeloupGoldbeterCircadianClock;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class LeloupGoldbeterCircadianClockSde implements Sde {

	public final double Na = 6.022 * 1e23;

	private LeloupGoldbeterCircadianClock ssaModel;

	private DoubleMatrix2D stoichiometryMatrix;
	private double omega;
	private double sqrtOmega;
	private double[] propensities;


	public LeloupGoldbeterCircadianClockSde(double volume, double inputOffset,
			double inputAmplitude, double inputFrequency) {
		this(volume, new SinusoidalFunction(inputOffset, inputAmplitude, inputFrequency));
	}

	public LeloupGoldbeterCircadianClockSde(double volume, ScalarTrajectoryFunction inputFunction) {
		omega = Na * volume;
		System.out.println("omega=" + omega);
		ssaModel = new LeloupGoldbeterCircadianClock(omega, inputFunction);
	    this.sqrtOmega = Math.sqrt(omega);
	    double[][] stoch = new double[ssaModel.getNumberOfReactions()][ssaModel.getNumberOfSpecies()];
		for (int i=0; i < ssaModel.getNumberOfReactions(); i++) {
			ssaModel.changeState(i, 0.0, stoch[i]);
		}
		stoichiometryMatrix = new DenseDoubleMatrix2D(stoch);
		propensities = new double[ssaModel.getNumberOfReactions()];
	}

	public LeloupGoldbeterCircadianClock getSSAModel() {
		return ssaModel;
	}

	public int getDriftDimension() {
		return ssaModel.getNumberOfSpecies();
	}

	public int getDiffusionDimension() {
		return ssaModel.getNumberOfReactions();
	}

	public double getOmega() {
		return omega;
	}

	public void computeDriftAndDiffusion(double t, DoubleMatrix1D X,
			DoubleMatrix1D F, DoubleMatrix2D G) throws Exception {
		double[] x = X.toArray();
		ssaModel.computePropensities(t, x, propensities);

		F.assign(0.0);
		G.assign(0.0);
		for (int k=0; k < getDiffusionDimension(); k++) {
			double propensity = propensities[k];
			double propensitySqrt = Math.sqrt(propensity);
			for (int i=0; i < getDriftDimension(); i++) {
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
