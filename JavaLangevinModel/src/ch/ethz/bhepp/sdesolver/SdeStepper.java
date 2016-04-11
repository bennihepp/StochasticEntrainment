package ch.ethz.bhepp.sdesolver;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public interface SdeStepper {

	Sde getSde();
	double getStepSize();

	double getTime();
	DoubleMatrix1D getState();
	double getState(int i);

	void setTime(double t);
	void setState(DoubleMatrix1D X);
	void setState(int i, double x);

	void step() throws Exception;
	void steps(int numOfSteps) throws Exception;

	void addStepHook(StepHook stepHook);

	interface StepHook {

		void call(double t, DoubleMatrix1D X, DoubleMatrix1D W, DoubleMatrix1D F, DoubleMatrix2D G);

	};

}
