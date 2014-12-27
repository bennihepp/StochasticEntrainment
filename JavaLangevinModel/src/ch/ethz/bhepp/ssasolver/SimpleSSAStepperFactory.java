package ch.ethz.bhepp.ssasolver;

import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

public class SimpleSSAStepperFactory implements SSAStepperFactory {

	private StochasticReactionNetworkModelFactory modelFactory;
	private RandomEngine rng;

	public SimpleSSAStepperFactory(StochasticReactionNetworkModelFactory modelFactory) {
		this.modelFactory = modelFactory;
	}

	public SimpleSSAStepperFactory(StochasticReactionNetworkModelFactory modelFactory, int seed) {
		this(modelFactory, new MersenneTwister(seed));
	}

	public SimpleSSAStepperFactory(StochasticReactionNetworkModelFactory modelFactory, RandomEngine rng) {
		this.modelFactory = modelFactory;
		this.rng = rng;
	}

	@Override
	public SSAStepper createStepper() {
		StochasticReactionNetworkModel model = modelFactory.createModel();
		return new SimpleSSAStepper(model, rng.nextInt());
	}

}
