package ch.ethz.bhepp.utils;
import cern.colt.matrix.DoubleMatrix1D;

public interface ScalarTrajectoryFunction {

	public double evaluate(double t, DoubleMatrix1D X);

	public double evaluate(double t, double[] X);

}
