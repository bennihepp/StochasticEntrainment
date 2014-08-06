package ch.ethz.bhepp.sdesolver;
import cern.colt.matrix.DoubleMatrix1D;

public interface ScalarTrajectoryFunction {

	public double evaluate(double t, DoubleMatrix1D X);

}
