package ch.ethz.bhepp.sdesolver;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class SdeSolution {

	private DoubleMatrix1D T;
	private DoubleMatrix2D X;

	public SdeSolution(DoubleMatrix1D T, DoubleMatrix2D X) {
		this.T = new DenseDoubleMatrix1D(T.toArray());
		this.X = new DenseDoubleMatrix2D(X.toArray());
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
