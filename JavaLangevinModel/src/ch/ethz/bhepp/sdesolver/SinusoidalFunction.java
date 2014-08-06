package ch.ethz.bhepp.sdesolver;
import cern.colt.matrix.DoubleMatrix1D;


public class SinusoidalFunction implements ScalarTrajectoryFunction {

	private double offset;
	private double amplitude;
	private double frequency;

	public SinusoidalFunction(double offset, double amplitude, double frequency) {
		this.offset = offset;
		this.amplitude = amplitude;
		this.frequency = frequency;
	}

	public double evaluate(double t, DoubleMatrix1D X) {
		return offset + amplitude * Math.sin(2 * Math.PI * t * frequency);
	}

}
