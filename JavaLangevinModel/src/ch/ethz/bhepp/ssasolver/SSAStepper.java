package ch.ethz.bhepp.ssasolver;
import cern.colt.matrix.DoubleMatrix1D;

public interface SSAStepper {

	StochasticReactionNetworkModel getModel();

	double getTime();
	DoubleMatrix1D getState();
	double getState(int i);

	void setTime(double t);
	void setState(DoubleMatrix1D X);
	void setState(int i, double x);

	void step() throws Exception;
	void steps(int numOfSteps) throws Exception;

}
