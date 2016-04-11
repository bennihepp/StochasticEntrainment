package ch.ethz.bhepp.ssasolver;

public interface StochasticReactionNetworkModel extends ReactionNetworkModel {

	double computePropensity(int reaction, double t, double[] x);

	void computePropensities(double t, double[] x, double[] propensities);

	double computePropensitiesAndSum(double t, double[] x, double[] propensities);

	void changeState(int reaction, double t, double[] x);

}
