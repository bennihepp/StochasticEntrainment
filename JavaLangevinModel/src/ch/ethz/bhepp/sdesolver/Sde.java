package ch.ethz.bhepp.sdesolver;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public interface Sde {

	int getDriftDimension();

	int getDiffusionDimension();
	
	/** Computes the drift F and diffusion G coefficients at a time t. */
	void computeDriftAndDiffusion(final double t, final DoubleMatrix1D X, 
			DoubleMatrix1D F, DoubleMatrix2D G) throws Exception;

}
