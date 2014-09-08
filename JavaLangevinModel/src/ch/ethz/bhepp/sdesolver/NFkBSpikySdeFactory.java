package ch.ethz.bhepp.sdesolver;

import ch.ethz.bhepp.sdesolver.models.NFkBLangevinSde;

public class NFkBSpikySdeFactory implements SdeFactory {

	private double volume;
	private ScalarTrajectoryFunction inputFunction;

	public NFkBSpikySdeFactory(double volume, ScalarTrajectoryFunction inputFunction) {
		this.volume = volume;
		this.inputFunction = inputFunction;
	}

	@Override
	public Sde createSde() {
		return new NFkBLangevinSde(volume, inputFunction);
	}

}
