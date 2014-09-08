package ch.ethz.bhepp.sdesolver;

import java.util.LinkedList;
import java.util.List;

import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;
import ch.ethz.bhepp.sdesolver.SdeStepper.StepHook;

public class EulerMaruyamaStepperFactory implements SdeStepperFactory {

	private SdeFactory sdeFactory;
	private double dt;
	private RandomEngine rng;
	private List<StepHook> stepHooks;

	public EulerMaruyamaStepperFactory(SdeFactory sdeFactory, double dt) {
		this(sdeFactory, dt, new MersenneTwister(new java.util.Date()));
	}

	public EulerMaruyamaStepperFactory(SdeFactory sdeFactory, double dt, int seed) {
		this(sdeFactory, dt, new MersenneTwister(seed));
	}

	public EulerMaruyamaStepperFactory(SdeFactory sdeFactory, double dt, RandomEngine rng) {
		this.sdeFactory = sdeFactory;
		this.dt = dt;
		this.rng = rng;
		stepHooks = new LinkedList<StepHook>();
	}

	public void addStepHook(StepHook stepHook) {
		stepHooks.add(stepHook);
	}

	public void clearStepHook() {
		stepHooks.clear();
	}

	@Override
	public SdeStepper createStepper() {
		Sde sde = sdeFactory.createSde();
		EulerMaruyamaStepper stepper = new EulerMaruyamaStepper(sde, dt, rng.nextInt());
		for (StepHook stepHook : stepHooks)
			stepper.addStepHook(stepHook);
		return stepper;
	}

}
