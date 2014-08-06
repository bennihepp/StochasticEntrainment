package ch.ethz.bhepp.sdesolver;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class MatrixHelper {

	public static DoubleMatrix1D createZeroDoubleMatrix1D(int size) {
		return new DenseDoubleMatrix1D(size);
	}

	public static DoubleMatrix1D createDoubleMatrix1D(double[] matrix) {
		return new DenseDoubleMatrix1D(matrix);
	}

	public static DoubleMatrix2D createZeroDoubleMatrix2D(int rows, int cols) {
		return new DenseDoubleMatrix2D(rows, cols);
	}

	public static DoubleMatrix2D createDoubleMatrix2D(double[][] matrix) {
		return new DenseDoubleMatrix2D(matrix);
	}

	public static double[] getMatrixData(DoubleMatrix1D matrix) {
		return matrix.toArray();
	}

	public static double[][] getMatrixData(DoubleMatrix2D matrix) {
		return matrix.toArray();
	}

}
