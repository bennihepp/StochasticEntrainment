package ch.ethz.bhepp.utils;
import cern.colt.matrix.DoubleMatrix1D;


public class SinusoidalFunction implements ScalarTrajectoryFunction {

	private double offset;
	private double amplitude;
	private double frequency;
	private double initial_phase;

	public SinusoidalFunction(double offset, double amplitude, double frequency) {
		this(offset, amplitude, frequency, 0.0);
	}

	public SinusoidalFunction(double offset, double amplitude, double frequency, double initial_phase) {
		this.offset = offset;
		this.amplitude = amplitude;
		this.frequency = frequency;
		this.initial_phase = initial_phase;
	}

	public double evaluate(double t) {
		return offset + amplitude * Math.sin(2 * Math.PI * t * frequency + initial_phase);
	}

	@Override
	public double evaluate(double t, DoubleMatrix1D X) {
		return evaluate(t);
	}

	@Override
	public double evaluate(double t, double[] X) {
		return evaluate(t);
	}

}
