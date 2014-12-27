package ch.ethz.bhepp.sdesolver.test;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.random.engine.MersenneTwister;
import ch.ethz.bhepp.sdesolver.EnsurePositiveStateHook;
import ch.ethz.bhepp.sdesolver.EulerMaruyamaStepper;
import ch.ethz.bhepp.sdesolver.SdeSolver;
import ch.ethz.bhepp.sdesolver.SdeStepper;
import ch.ethz.bhepp.sdesolver.models.CircadianClockDrosophilaSde;
import ch.ethz.bhepp.utils.FiniteTimeSolution;
import ch.ethz.bhepp.utils.MatrixHelper;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class CircadianClockDrosophilaSolverTest {

	public static void main(String[] args) throws Exception {
		double t0 = 0.0;
		double tf = 1000.0;
		double volume = 1e-15;
		double dt = 0.0001;

		ScalarTrajectoryFunction inputFunction = new SinusoidalFunction(1.0, 0.14, 1 / 3.0);
		CircadianClockDrosophilaSde model = new CircadianClockDrosophilaSde(volume, inputFunction);

		SdeStepper stepper = new EulerMaruyamaStepper(model, dt, new MersenneTwister(0));
		stepper.addStepHook(new EnsurePositiveStateHook());
//		stepper.addStepHook(new EnsureBoundedStateHook(0, 0, 1.0));
		SdeSolver solver = new SdeSolver(stepper);

		double[] x0 = { 0.0, 0.0, 0.0, 0.0, 0.0 };
		DoubleMatrix1D X0 = MatrixHelper.createDoubleMatrix1D(x0);
//		DoubleMatrix1D X0 = MatrixHelper.createZeroDoubleMatrix1D(model.getDriftDimension());
		FiniteTimeSolution solution = solver.solve(t0, X0, tf);

		DoubleMatrix1D T = solution.getT();
		DoubleMatrix2D X = solution.getX();

		System.out.println("T.size()=" + T.size());
		System.out.println("X.rows()=" + X.rows() + " X.columns()=" + X.columns());

		System.out.println(X.get(X.rows() - 1, 0));
	}

}
