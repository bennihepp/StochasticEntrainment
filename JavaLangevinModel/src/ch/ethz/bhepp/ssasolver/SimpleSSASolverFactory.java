package ch.ethz.bhepp.ssasolver;

public class SimpleSSASolverFactory implements SSASolverFactory {

	private SSAStepperFactory stepperFactory;
	private int recordIndex;

	public SimpleSSASolverFactory(SSAStepperFactory stepperFactory) {
		this(stepperFactory, -1);
	}

	public SimpleSSASolverFactory(SSAStepperFactory stepperFactory, int recordIndex) {
		this.stepperFactory = stepperFactory;
		this.recordIndex = recordIndex;
	}

	@Override
	public SSASolver createSolver() {
		SSAStepper stepper = stepperFactory.createStepper();
		return new SSASolver(stepper, recordIndex);
	}

}
