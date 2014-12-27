package ch.ethz.bhepp.ssasolver.test;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.random.engine.MersenneTwister;
import ch.ethz.bhepp.ssasolver.SSASolver;
import ch.ethz.bhepp.ssasolver.SSAStepper;
import ch.ethz.bhepp.ssasolver.SimpleSSAStepper;
import ch.ethz.bhepp.ssasolver.models.LeloupGoldbeterCircadianClock;
import ch.ethz.bhepp.utils.FiniteTimeSolution;
import ch.ethz.bhepp.utils.MatrixHelper;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class LeloupGoldbeterCircadianClockTest {

	public static void main(String[] args) throws Exception {
		double t0 = 0.0;
		double tf = 100.0;
		double recordStep = tf / 1000.0;
		double omega = 600;

		ScalarTrajectoryFunction inputFunction = new SinusoidalFunction(0.0, 0.0, 1.0);
		LeloupGoldbeterCircadianClock model = new LeloupGoldbeterCircadianClock(omega, inputFunction);

		SSAStepper stepper = new SimpleSSAStepper(model, new MersenneTwister(0));
		SSASolver solver = new SSASolver(stepper);

		double[] x0 = new double[model.getNumberOfSpecies()];
		for (int i=0; i < x0.length; i++)
			x0[i] = 0.0;
		DoubleMatrix1D X0 = MatrixHelper.createDoubleMatrix1D(x0);
//		DoubleMatrix1D X0 = MatrixHelper.createZeroDoubleMatrix1D(model.getDriftDimension());
		FiniteTimeSolution solution = solver.solve(t0, X0, tf, recordStep);

		DoubleMatrix1D T = solution.getT();
		DoubleMatrix2D X = solution.getX();

		System.out.println("T.size()=" + T.size());
		System.out.println("X.rows()=" + X.rows() + " X.columns()=" + X.columns());

		for (int i=0; i < X.rows(); i++)
			System.out.println(X.get(i, 0));
	}

}
