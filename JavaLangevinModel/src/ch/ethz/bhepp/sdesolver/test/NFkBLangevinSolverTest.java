package ch.ethz.bhepp.sdesolver.test;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.random.engine.MersenneTwister;
import ch.ethz.bhepp.sdesolver.EnsurePositiveStateHook;
import ch.ethz.bhepp.sdesolver.EulerMaruyamaStepper;
import ch.ethz.bhepp.sdesolver.SdeSolver;
import ch.ethz.bhepp.sdesolver.SdeStepper;
import ch.ethz.bhepp.sdesolver.models.NFkBLangevinSde;
import ch.ethz.bhepp.utils.FiniteTimeSolution;
import ch.ethz.bhepp.utils.MatrixHelper;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;
import ch.ethz.bhepp.utils.SinusoidalFunction;

public class NFkBLangevinSolverTest {

	public static void main(String[] args) throws Exception {
		double t0 = 0.0;
		double tf = 20000.0;
		double volume = 1e-13;
		double dt = 0.01;

		ScalarTrajectoryFunction inputFunction = new SinusoidalFunction(0.1, 0.05, 1 / 120.0);
		NFkBLangevinSde model = new NFkBLangevinSde(volume, inputFunction);
		model.setN_tot(1.2203);
		model.setI_tot(1.0070);
		model.setK_Nin(3.6552);
		model.setK_I(0.0289);
		model.setK_lin(0.0236);
		model.setK_N(0.0265);
		model.setK_t(0.7456);
		model.setGamma_m(0.0163);
		model.setK_tl(0.3302);
		model.setAlpha(0.8811);
		model.setK_a(0.2354);
		model.setK_i(0.1762);
		model.setK_p(0.0360);
		model.setK_A20(0.0018);
		model.setA_20(0.0035);

		SdeStepper stepper = new EulerMaruyamaStepper(model, dt, new MersenneTwister(0));
		stepper.addStepHook(new EnsurePositiveStateHook());
		SdeSolver solver = new SdeSolver(stepper);

		double[] x0 = { 0.3300, 1.0000, 2.7500, 0.1400, 0.9400 };
		DoubleMatrix1D X0 = MatrixHelper.createDoubleMatrix1D(x0);
//		DoubleMatrix1D X0 = MatrixHelper.createZeroDoubleMatrix1D(model.getDriftDimension());
		FiniteTimeSolution solution = solver.solve(t0, X0, tf);

		DoubleMatrix1D T = solution.getT();
		DoubleMatrix2D X = solution.getX();

		System.out.println("T.size()=" + T.size());
		System.out.println("X.rows()=" + X.rows() + " X.columns()=" + X.columns());
	}

}
