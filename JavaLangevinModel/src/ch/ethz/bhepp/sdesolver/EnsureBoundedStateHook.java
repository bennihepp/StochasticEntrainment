package ch.ethz.bhepp.sdesolver;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.SdeStepper.StepHook;

public class EnsureBoundedStateHook implements StepHook {

	private int index;
	private double lowerBound;
	private double upperBound;

	public EnsureBoundedStateHook(int index, double lowerBound, double upperBound) {
		assert(lowerBound <= upperBound);
		this.index = index;
		this.lowerBound = lowerBound;
		this.upperBound = upperBound;
	}

	public void call(double t, DoubleMatrix1D X, DoubleMatrix1D W, DoubleMatrix1D F, DoubleMatrix2D G) {
		double value = X.get(index);
		if (value < lowerBound)
			X.set(index, lowerBound);
		else if (value > upperBound)
			X.set(index, upperBound);
	}

}
