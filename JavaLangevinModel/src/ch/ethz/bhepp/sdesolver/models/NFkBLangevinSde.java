package ch.ethz.bhepp.sdesolver.models;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.ethz.bhepp.sdesolver.ScalarTrajectoryFunction;
import ch.ethz.bhepp.sdesolver.Sde;
import ch.ethz.bhepp.sdesolver.SinusoidalFunction;

public class NFkBLangevinSde implements Sde {

	public final int NUM_OF_SPECIES = 5;
	public final int NUM_OF_REACTIONS = 9;

	public final double Na = 6.022 * 1e23;

	private double N_tot = 1.0;
	private double I_tot = 2.0;

	private double k_Nin = 5.4;
	private double K_I = 0.035;
	private double k_lin = 0.018;
	private double K_N = 0.029;
	private double k_t = 1.03;
	private double gamma_m = 0.017;
	private double k_tl = 0.24;
	private double alpha = 1.05;
	private double k_a = 0.24;
	private double k_i = 0.18;
	private double k_p = 0.036;
	private double k_A20 = 0.0018;
	private double A_20 = 0.0028;

	private DoubleMatrix2D stoichiometryMatrix;
	private DoubleMatrix1D scaledPropensities;
	private double omega;
	private double sqrtOmega;
	private ScalarTrajectoryFunction inputFunction;


	public NFkBLangevinSde(double volume, double inputOffset, double inputAmplitude, double inputFrequency) {
		this(volume, new SinusoidalFunction(inputOffset, inputAmplitude, inputFrequency));
	}

	public NFkBLangevinSde(double volume, ScalarTrajectoryFunction inputFunction) {
		omega = 1e-6 * Na * volume;
	    this.sqrtOmega = Math.sqrt(omega);
		this.inputFunction = inputFunction;
		stoichiometryMatrix = new DenseDoubleMatrix2D(NUM_OF_REACTIONS, NUM_OF_SPECIES);
		stoichiometryMatrix.set(0, 0, +1);
		stoichiometryMatrix.set(1, 0, -1);
		stoichiometryMatrix.set(2, 1, +1);
		stoichiometryMatrix.set(3, 1, -1);
		stoichiometryMatrix.set(4, 2, +1);
		stoichiometryMatrix.set(5, 2, -1);
		stoichiometryMatrix.set(6, 3, +1);
		stoichiometryMatrix.set(7, 3, -1);
		stoichiometryMatrix.set(7, 4, +1);
		stoichiometryMatrix.set(8, 4, -1);
		scaledPropensities = new DenseDoubleMatrix1D(NUM_OF_REACTIONS);
	}

	public int getDriftDimension() {
		return NUM_OF_SPECIES;
	}

	public int getDiffusionDimension() {
		return NUM_OF_REACTIONS;
	}

	public double getOmega() {
		return omega;
	}

	public void computeDriftAndDiffusion(double t, DoubleMatrix1D X,
			DoubleMatrix1D F, DoubleMatrix2D G) throws Exception {

		double TNF = inputFunction.evaluate(t, X);
		F.set(0, k_Nin * (N_tot - X.get(0)) * K_I/(K_I + X.get(2)) - k_lin*X.get(2)*X.get(0)/(K_N + X.get(0)));
		F.set(1, k_t*X.get(0)*X.get(0) - gamma_m*X.get(1));
		F.set(2, k_tl*X.get(1) - alpha*X.get(3)*(N_tot - X.get(0))*X.get(2)/(K_I + X.get(2)));
		F.set(3, k_a*TNF*(I_tot - X.get(3) - X.get(4)) - k_i*X.get(3));
		F.set(4, k_i*X.get(3) - k_p*X.get(4)*k_A20/(k_A20 + A_20*TNF));

        scaledPropensities.set(0, k_Nin * (N_tot - X.get(0)) * K_I / (K_I + X.get(2)));
        scaledPropensities.set(1, k_lin * X.get(2) * X.get(0) / (K_N + X.get(0)));
        scaledPropensities.set(2, k_t * X.get(0) * X.get(0));
        scaledPropensities.set(3, gamma_m * X.get(1));
        scaledPropensities.set(4, k_tl * X.get(1));
        scaledPropensities.set(5, alpha * X.get(3) * (N_tot - X.get(0)) * X.get(2) / (K_I + X.get(2)));
        scaledPropensities.set(6, k_a * TNF * (I_tot - X.get(3) - X.get(4)));
        scaledPropensities.set(7, k_i * X.get(3));
        scaledPropensities.set(8, k_p * X.get(4) * k_A20 / (k_A20 + A_20 * TNF));
//        for (int k=0; k < scaledPropensities.size(); k++)
//        	if (scaledPropensities.get(k) < 0.0)
//        		scaledPropensities.set(k, 0.0);

		G.assign(0.0);
		for (int k=0; k < NUM_OF_REACTIONS; k++) {
			double propensitySqrt = Math.sqrt(scaledPropensities.get(k));
			for (int i=0; i < NUM_OF_SPECIES; i++) {
				G.set(i, k, G.get(i, k) + propensitySqrt * stoichiometryMatrix.get(k, i) / sqrtOmega);
			}
		}

	}

	public double getN_tot() {
		return N_tot;
	}

	public void setN_tot(double n_tot) {
		this.N_tot = n_tot;
	}

	public double getI_tot() {
		return I_tot;
	}

	public void setI_tot(double i_tot) {
		I_tot = i_tot;
	}

	public double getK_Nin() {
		return k_Nin;
	}

	public void setK_Nin(double k_Nin) {
		this.k_Nin = k_Nin;
	}

	public double getK_I() {
		return K_I;
	}

	public void setK_I(double k_I) {
		K_I = k_I;
	}

	public double getK_lin() {
		return k_lin;
	}

	public void setK_lin(double k_lin) {
		this.k_lin = k_lin;
	}

	public double getK_N() {
		return K_N;
	}

	public void setK_N(double k_N) {
		K_N = k_N;
	}

	public double getK_t() {
		return k_t;
	}

	public void setK_t(double k_t) {
		this.k_t = k_t;
	}

	public double getGamma_m() {
		return gamma_m;
	}

	public void setGamma_m(double gamma_m) {
		this.gamma_m = gamma_m;
	}

	public double getK_tl() {
		return k_tl;
	}

	public void setK_tl(double k_tl) {
		this.k_tl = k_tl;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public double getK_a() {
		return k_a;
	}

	public void setK_a(double k_a) {
		this.k_a = k_a;
	}

	public double getK_i() {
		return k_i;
	}

	public void setK_i(double k_i) {
		this.k_i = k_i;
	}

	public double getK_p() {
		return k_p;
	}

	public void setK_p(double k_p) {
		this.k_p = k_p;
	}

	public double getK_A20() {
		return k_A20;
	}

	public void setK_A20(double k_A20) {
		this.k_A20 = k_A20;
	}

	public double getA_20() {
		return A_20;
	}

	public void setA_20(double a_20) {
		A_20 = a_20;
	}

}
