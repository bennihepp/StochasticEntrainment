package ch.ethz.bhepp.sdesolver;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.ethz.bhepp.utils.FiniteTimeSolution;

public class SdeSolver {

	private SdeStepper stepper;

	public SdeSolver(SdeStepper stepper) {
		this.stepper = stepper;
	}

	public FiniteTimeSolution solve(double t0, DoubleMatrix1D X0, double tf) throws Exception {
		return solve(t0, X0, tf, 1.0);
	}

	public FiniteTimeSolution solve(double t0, DoubleMatrix1D X0, double tf, double recordStep) throws Exception {
		int numOfStepsToRecord = (int)Math.round(Math.ceil((tf - t0) / recordStep) + 1);
		stepper.setTime(t0);
		stepper.setState(X0);
		DoubleMatrix1D T = new DenseDoubleMatrix1D(numOfStepsToRecord);
		DoubleMatrix2D X = new DenseDoubleMatrix2D(numOfStepsToRecord, X0.size());
		int k = 0;
		record(0, T, X);
		k++;
		while (stepper.getTime() <= tf || k < numOfStepsToRecord) {
			stepper.step();
			if (stepper.getTime() >= k * recordStep) {
				record(k, T, X);
				T.set(k, stepper.getTime());
				for (int l=0; l < X0.size(); l++) {
					X.set(k, l, stepper.getState(l));
				}
				k++;
			}
		}
		return new FiniteTimeSolution(T, X);
	}

	private void record(int k, DoubleMatrix1D T, DoubleMatrix2D X) {
		T.set(k, stepper.getTime());
		for (int l=0; l < X.columns(); l++) {
			X.set(k, l, stepper.getState(l));
		}
	}

}
