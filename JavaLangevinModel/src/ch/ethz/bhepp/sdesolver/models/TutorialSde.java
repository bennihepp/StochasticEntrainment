package ch.ethz.bhepp.sdesolver.models;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.Sde;

public class TutorialSde implements Sde {

	/** The number of species. */
	private final int NUM_OF_SPECIES = 2;
	
	/** The number of reactions. */
	private final int NUM_OF_REACTIONS = 4;

	private final double TRANSCRIPTION_RATE = 1.0;
	private final double MRNA_DEGRADATION_RATE = 0.1;
	private final double TRANSLATION_RATE = 0.1;
	private final double PROTEIN_DEGRADATION_RATE = 0.01;

	private double sqrtOmega = 1;

	@Override
	public int getDriftDimension() {
		return NUM_OF_SPECIES;
	}

	@Override
	public int getDiffusionDimension() {
		return NUM_OF_REACTIONS;
	}

	@Override
	public void computeDriftAndDiffusion(double t, DoubleMatrix1D X,
			DoubleMatrix1D F, DoubleMatrix2D G) throws Exception {
		double mRNA = X.get(0);
		double protein = X.get(1);
		double prop1 = TRANSCRIPTION_RATE;
		double prop2 = mRNA * MRNA_DEGRADATION_RATE;
		double prop3 = mRNA * TRANSLATION_RATE;
		double prop4 = protein * PROTEIN_DEGRADATION_RATE;
		F.set(0, prop1 - prop2);
		F.set(1, prop3 - prop4);

		G.assign(0.0);
		double propSqrt1 = Math.sqrt(prop1);
		double propSqrt2 = Math.sqrt(prop2);
		double propSqrt3 = Math.sqrt(prop3);
		double propSqrt4 = Math.sqrt(prop4);
		G.set(0, 0, propSqrt1 / sqrtOmega);
		G.set(0, 1, propSqrt2 / sqrtOmega);
		G.set(1, 2, propSqrt3 / sqrtOmega);
		G.set(1, 3, propSqrt4 / sqrtOmega);
	}

}
