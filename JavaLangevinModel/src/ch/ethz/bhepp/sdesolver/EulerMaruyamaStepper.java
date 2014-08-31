package ch.ethz.bhepp.sdesolver;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.random.Normal;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;

public class EulerMaruyamaStepper implements SdeStepper {

	private Sde sde;
	private double dt;
	private double sqrtDt;
	private Normal normalDistribution;

	private double _t;
	private DoubleMatrix1D _X;
	private DoubleMatrix1D _W;
	private DoubleMatrix1D _F;
	private DoubleMatrix2D _G;
	private List<StepHook> stepHooks;

	public EulerMaruyamaStepper(Sde sde, double dt) {
		this(sde, dt, new MersenneTwister(new java.util.Date()));
	}

	public EulerMaruyamaStepper(Sde sde, double dt, int seed) {
		this(sde, dt, new MersenneTwister(seed));
	}

	public EulerMaruyamaStepper(Sde sde, double dt, RandomEngine rng) {
		this.sde = sde;
		this.dt = dt;
		sqrtDt = Math.sqrt(dt);
		normalDistribution = new Normal(0.0, 1.0, rng);
		_X = new DenseDoubleMatrix1D(sde.getDriftDimension());
		_W = new DenseDoubleMatrix1D(sde.getDiffusionDimension());
		_F = new DenseDoubleMatrix1D(sde.getDriftDimension());
		_G = new DenseDoubleMatrix2D(sde.getDriftDimension(), sde.getDiffusionDimension());
		stepHooks = new LinkedList<StepHook>();
	}

	public void clearStepHook() {
		stepHooks.clear();
	}

	public void addStepHook(StepHook stepHook) {
		this.stepHooks.add(stepHook);
	}

	public Sde getSde() {
		return sde;
	}

	public double getStepSize() {
		return dt;
	}

	public double getTime() {
		return _t;
	}

	public DoubleMatrix1D getState() {
		return new DenseDoubleMatrix1D(_X.toArray());
	}

	public double getState(int i) {
		return _X.get(i);
	}

	public void setTime(double t) {
		_t = t;
	}

	public void setState(DoubleMatrix1D X) {
		_X.assign(X);
	}

	public void setState(int i, double x) {
		_X.set(i, x);
	}

	public void step() throws Exception {
		sde.computeDriftAndDiffusion(_t, _X, _F, _G);
		for (int k=0; k < _W.size(); k++) {
			_W.set(k, normalDistribution.nextDouble());
		}
		for (int i=0; i < _X.size(); i++) {
			double x = _X.get(i) + _F.get(i) * dt;
			double noise = 0.0;
			for (int k=0; k < _W.size(); k++) {
				noise += _G.get(i, k) * _W.get(k);
			}
			double newValue = x + noise * sqrtDt;
			_X.set(i, newValue);
		}
		Iterator<StepHook> it = stepHooks.iterator();
		while (it.hasNext()) {
			StepHook stepHook = it.next();
			stepHook.call(_t, _X, _W, _F, _G);
		}
		_t += dt;
	}

	public void steps(int numOfSteps) throws Exception {
		for (int i=0; i < numOfSteps; i++)
			step();
	}

}
