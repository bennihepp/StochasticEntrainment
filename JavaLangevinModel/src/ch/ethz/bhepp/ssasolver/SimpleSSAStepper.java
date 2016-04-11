package ch.ethz.bhepp.ssasolver;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.jet.random.Empirical;
import cern.jet.random.EmpiricalWalker;
import cern.jet.random.Exponential;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;


public class SimpleSSAStepper implements SSAStepper {

	private StochasticReactionNetworkModel model;
//	private Uniform discreteUniform;
	private EmpiricalWalker empiricalWalker;
	private Exponential exponential;
	private double t;
	private double[] X;
	private double[] propVec;

	public SimpleSSAStepper(StochasticReactionNetworkModel model) {
		this(model, new MersenneTwister(new java.util.Date()));
	}

	public SimpleSSAStepper(StochasticReactionNetworkModel model, int seed) {
		this(model, new MersenneTwister(seed));
	}

	public SimpleSSAStepper(StochasticReactionNetworkModel model, RandomEngine rng) {
		this.model = model;
		this.t = 0.0;
		this.X = new double[model.getNumberOfSpecies()];
		this.propVec = new double[model.getNumberOfReactions()];
		propVec[0] = 1.0;
		this.empiricalWalker = new EmpiricalWalker(propVec, Empirical.NO_INTERPOLATION, rng);
//		this.discreteUniform = new Uniform(0, model.getNumberOfReactions() - 1, rng);
		this.exponential = new Exponential(1.0, rng);
	}

//	public double simulate(StochasticReactionNetworkModel model, double t0, double[] x0, double t1, double[] x1) {
//		assert(x0.length == x1.length);
//		assert(x0.length == model.getNumberOfSpecies());
//
//		final int STOCHASTIC_RECORD_INTERVAL = 1;
//
//		double[] x = new double[x0.length];
//		for (int i=0; i < x0.length; i++)
//			x[i] = x0[i];
//		double t = t0;
//		double[] propVec = new double[model.getNumberOfReactions()];
//    	for (TrajectoryRecorder handler : trajectoryRecorders) {
//    		handler.beginRecording(t0, x0, t1);
//    	}
//    	long reactionCounter = 0;
//    	double[] reactionCounterArray = new double[model.getNumberOfReactions()];
//    	double msgDt = (t1 - t0) / 20.0;
//    	double nextMsgT = t0 + msgDt;
//		final long startTime = System.currentTimeMillis();
//
//		int j = 0;
//
//		while (true) {
//
//			if (showProgress)
//				while (t > nextMsgT) {
//					System.out.println("Progress: " + (100 * (nextMsgT - t0)/(t1 - t0)) + "%");
//					nextMsgT += msgDt;
//				}
//
//	        double propSum = model.computePropensitiesAndSum(t, x, propVec);
//
//	        // Sanity check for negative propensities or no more reactions occuring (zero propensities)
//	        // FIXME: NaN check necessary?
//	        if (propSum < 0 || Double.isNaN(propSum)) {
//	        	throw new RuntimeException("Negative propensities are not allowed to occur!");
//	        }
//
////	        if (propSum == 0) {
////	        	t = t1;
////	        	break;
////	        }
//
//	        double tau = rdg.nextExponential(1 / propSum);
//	        t = t + tau;
//
//	        // Stop if we reached the end-timepoint
//	        if (t >= t1)
//	        	break;
//
//	        // Determine which reaction fired and update state
//	        int reaction = DiscreteProbabilityDistribution.sample(rdg, propVec, propSum);
//    		model.changeState(reaction, t, x);
//
//        	reactionCounter++;
//        	reactionCounterArray[reaction]++;
//
//        	j++;
//        	if (j >= STOCHASTIC_RECORD_INTERVAL) {
//	        	for (TrajectoryRecorder handler : trajectoryRecorders)
//	        		handler.record(t, x);
//        		j = 0;
//        	}
//
//		}
//
//		final long endTime = System.currentTimeMillis();
//		if (printMessages) {
//			System.out.println("Execution time: " + (endTime - startTime));
//			System.out.println("Total of " + reactionCounter + " reactions performed");
//			Utilities.printArray("Total reaction counts", reactionCounterArray);
//			for (int r=0; r < reactionCounterArray.length; r++)
//				reactionCounterArray[r] = 100.0 * reactionCounterArray[r] / reactionCounter;
//			Utilities.printArray("Relative reaction counts", reactionCounterArray, "%.2f");
//		}
//		for (int i=0; i < x1.length; i++)
//			x1[i] = x[i];
//    	for (TrajectoryRecorder handler : trajectoryRecorders)
//    		handler.endRecording(x1);
//		return t;
//	}

	@Override
	public StochasticReactionNetworkModel getModel() {
		return model;
	}

	@Override
	public double getTime() {
		return t;
	}

	@Override
	public DoubleMatrix1D getState() {
		return new DenseDoubleMatrix1D(X);
	}

	@Override
	public double getState(int i) {
		return X[i];
	}

	@Override
	public void setTime(double t) {
		this.t = t;
	}

	@Override
	public void setState(DoubleMatrix1D X) {
		for (int i=0; i < this.X.length; i++)
			this.X[i] = X.get(i);
	}

	@Override
	public void setState(int i, double x) {
		this.X[i] = x;
	}

	@Override
	public void step() throws Exception {
        double propSum = model.computePropensitiesAndSum(t, X, propVec);

        // Sanity check for negative propensities or no more reactions occuring (zero propensities)
        // FIXME: NaN check necessary?
        if (propSum < 0 || Double.isNaN(propSum)) {
        	throw new RuntimeException("Negative propensities are not allowed to occur!");
        }

        if (propSum == 0) {
        	t = Double.POSITIVE_INFINITY;
        	return;
        }

        double tau = exponential.nextDouble(propSum);
        t = t + tau;

        // Determine which reaction fired and update state
        empiricalWalker.setState2(propVec);
        int reaction = empiricalWalker.nextInt();
//        int reaction = discreteUniform.nextIntFromTo(0, model.getNumberOfReactions() - 1);
		model.changeState(reaction, t, X);
	}

	@Override
	public void steps(int numOfSteps) throws Exception {
		for (int i=0; i < numOfSteps; i++)
			step();
	}

}
