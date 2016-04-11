package ch.ethz.bhepp.ssasolver.models;

import ch.ethz.bhepp.ssasolver.StochasticReactionNetworkModel;
import ch.ethz.bhepp.ssasolver.StochasticReactionNetworkModelFactory;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;

public class LeloupGoldbeterCircadianClockFactory implements StochasticReactionNetworkModelFactory {

	private double omega;
	private ScalarTrajectoryFunction inputFunction;

	public LeloupGoldbeterCircadianClockFactory(double omega, ScalarTrajectoryFunction inputFunction) {
		this.omega = omega;
		this.inputFunction = inputFunction;
	}

	@Override
	public StochasticReactionNetworkModel createModel() {
		return new LeloupGoldbeterCircadianClock(omega, inputFunction);
	}

}
