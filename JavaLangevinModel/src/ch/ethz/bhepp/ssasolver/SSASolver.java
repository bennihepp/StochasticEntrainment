package ch.ethz.bhepp.ssasolver;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.ethz.bhepp.utils.FiniteTimeSolution;

public class SSASolver {

	private SSAStepper stepper;
	private int recordIndex;

	public SSASolver(SSAStepper stepper) {
		this(stepper, -1);
	}

	public SSASolver(SSAStepper stepper, int recordIndex) {
		this.stepper = stepper;
		this.recordIndex = recordIndex;
	}

	public FiniteTimeSolution solve(double t0, DoubleMatrix1D X0, double tf) throws Exception {
		return solve(t0, X0, tf, 1.0);
	}

	public FiniteTimeSolution solve(double t0, DoubleMatrix1D X0, double tf, double recordStep) throws Exception {
		int numOfStepsToRecord = (int)Math.round(Math.ceil((tf - t0) / recordStep) + 1);
		stepper.setTime(t0);
		stepper.setState(X0);
		DoubleMatrix1D T = new DenseDoubleMatrix1D(numOfStepsToRecord);
		int sizeX = recordIndex >=0 ? 1 : X0.size();
		DoubleMatrix2D X = new DenseDoubleMatrix2D(numOfStepsToRecord, sizeX);
		int k = 0;
		record(k, t0, X0, T, X);
		k++;
		while (stepper.getTime() <= tf || k < numOfStepsToRecord) {
			stepper.step();
//			System.out.println(stepper.getState(0));
			while (stepper.getTime() > k * recordStep && k < X.rows()) {
				double t = k * recordStep;
				DoubleMatrix1D x = stepper.getState();
				record(k, t, x, T, X);
//				T.set(k, t);
//				for (int l=0; l < X0.size(); l++) {
//					X.set(k, l, X.get(k-1, l));
//				}
				k++;
			}
		}
		return new FiniteTimeSolution(T, X);
	}

	private void record(int k, double t, DoubleMatrix1D x, DoubleMatrix1D T, DoubleMatrix2D X) {
		T.set(k, t);
		if (recordIndex >=0) {
			X.set(k, 0, x.get(recordIndex));
		} else {
			for (int l=0; l < x.size(); l++) {
				X.set(k, l, x.get(l));
			}
		}
	}

}
