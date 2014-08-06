package ch.ethz.bhepp.sdesolver.test;

import cern.jet.random.engine.MersenneTwister;
import ch.ethz.bhepp.sdesolver.EnsurePositiveStateHook;
import ch.ethz.bhepp.sdesolver.EulerMaruyamaStepper;
import ch.ethz.bhepp.sdesolver.ScalarTrajectoryFunction;
import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.sdesolver.SdeStepper;
import ch.ethz.bhepp.sdesolver.SinusoidalFunction;
import ch.ethz.bhepp.sdesolver.models.NFkBLangevinSde;

public class NFkBLangevinStepperTest {

	public static void main(String[] args) throws Exception {
		double tf = 20000.0;
		double volume = 1e-13;
		double dt = 0.01;
		long numOfSteps = Math.round(tf / dt);

		ScalarTrajectoryFunction inputFunction = new SinusoidalFunction(0.1, 0.05, 1 / 120.0);
		Sde model = new NFkBLangevinSde(volume, inputFunction);
		SdeStepper stepper = new EulerMaruyamaStepper(model, dt, new MersenneTwister(0));
		stepper.addStepHook(new EnsurePositiveStateHook());

		System.out.println(stepper.getTime());
		System.out.println(stepper.getState());

		for (long i=0; i < numOfSteps; i++) {
			stepper.step();
		}

		System.out.println(stepper.getTime());
		System.out.println(stepper.getState());

	}

}
