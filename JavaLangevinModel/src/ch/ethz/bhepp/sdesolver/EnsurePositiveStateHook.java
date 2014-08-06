package ch.ethz.bhepp.sdesolver;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.SdeStepper.StepHook;

public class EnsurePositiveStateHook implements StepHook {

	public void call(double t, DoubleMatrix1D X, DoubleMatrix1D W, DoubleMatrix1D F, DoubleMatrix2D G) {
		for (int i=0; i < X.size(); i++) {
			if (X.get(i) < 0.0)
				X.set(i, 0.0);
		}
	}

}
