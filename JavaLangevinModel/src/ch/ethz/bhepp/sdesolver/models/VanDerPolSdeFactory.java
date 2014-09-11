package ch.ethz.bhepp.sdesolver.models;

import ch.ethz.bhepp.sdesolver.ScalarTrajectoryFunction;
import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.sdesolver.SdeFactory;

public class VanDerPolSdeFactory implements SdeFactory {

	private double volume;
	private ScalarTrajectoryFunction inputFunction;

	public VanDerPolSdeFactory(double volume, ScalarTrajectoryFunction inputFunction) {
		this.volume = volume;
		this.inputFunction = inputFunction;
	}

	@Override
	public Sde createSde() {
		return new VanDerPolSde(volume, inputFunction);
	}

}
