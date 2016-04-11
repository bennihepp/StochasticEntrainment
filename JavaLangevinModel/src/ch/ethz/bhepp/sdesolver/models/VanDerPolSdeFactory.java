package ch.ethz.bhepp.sdesolver.models;

import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.sdesolver.SdeFactory;
import ch.ethz.bhepp.utils.ScalarTrajectoryFunction;

public class VanDerPolSdeFactory implements SdeFactory {

	private double volume;
	private ScalarTrajectoryFunction inputFunction;
	private boolean setParameters = false;
	private double B;
	private double d;

	public VanDerPolSdeFactory(double volume, ScalarTrajectoryFunction inputFunction) {
		this.volume = volume;
		this.inputFunction = inputFunction;
	}

	public VanDerPolSdeFactory(double volume, ScalarTrajectoryFunction inputFunction, double B, double d) {
		this.volume = volume;
		this.inputFunction = inputFunction;
		this.setParameters = true;
		this.B = B;
		this.d = d;
	}

	@Override
	public Sde createSde() {
		VanDerPolSde sde = new VanDerPolSde(volume, inputFunction);
		if (this.setParameters) {
			sde.setParameterB(this.B);
			sde.setParameterD(this.d);
		}
		return sde;
	}

}
