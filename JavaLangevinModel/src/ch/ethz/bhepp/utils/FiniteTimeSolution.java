package ch.ethz.bhepp.utils;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class FiniteTimeSolution {

	private DoubleMatrix1D T;
	private DoubleMatrix2D X;

	public FiniteTimeSolution(DoubleMatrix1D T, DoubleMatrix2D X) {
		this.T = T;
		this.X = X;
	}

	public DoubleMatrix1D getT() {
		return new DenseDoubleMatrix1D(T.toArray());
	}

	public DoubleMatrix2D getX() {
		return new DenseDoubleMatrix2D(X.toArray());
	}

	public double[] getTArray() {
		return T.toArray();
	}

	public double[][] getXArray() {
		return X.toArray();
	}

}
