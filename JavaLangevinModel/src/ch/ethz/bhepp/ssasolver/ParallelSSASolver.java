package ch.ethz.bhepp.ssasolver;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import ch.ethz.bhepp.utils.FiniteTimeSolution;

public class ParallelSSASolver {

	private static final int DEFAULT_NUM_OF_THREADS = 4;
	private int numOfThreads;
	private ExecutorService executor;
	private SSASolverFactory solverFactory;

	public ParallelSSASolver(SSASolverFactory solverFactory) {
		this(solverFactory, 0);
	}

	public ParallelSSASolver(SSASolverFactory solverFactory, int numOfThreads) {
		this.solverFactory = solverFactory;
		if (numOfThreads == 0) {
			String numOfThreadsStr = System.getenv("NUM_OF_THREADS");
			boolean foundNumOfThreads = false;
			if (numOfThreadsStr != null) {
				try {
					numOfThreads = Integer.parseInt(numOfThreadsStr);
					foundNumOfThreads = true;
				} catch (NumberFormatException e) {
					System.out.println("WARNING: Couldn't parse environment variable NUM_OF_THREADS");
				}
			}
			if (!foundNumOfThreads) {
				numOfThreads = DEFAULT_NUM_OF_THREADS;
				System.out.println("WARNING: Using default number of threads: " + numOfThreads);
			}
		}
		this.numOfThreads = numOfThreads;
	}

	public void startPool() {
		if (executor != null) {
			shutdownPool();
		}
		executor = Executors.newFixedThreadPool(numOfThreads);
	}

	public void shutdownPool() {
		executor.shutdown();
		executor = null;
	}

	private void ensurePoolStarted() {
		if (executor == null)
			startPool();
	}

	public List<FiniteTimeSolution> solve(int samples, double t0, DoubleMatrix2D X0, double tf) throws Exception {
		return solve(samples, t0, X0, tf, 1.0, false);
	}

	public List<FiniteTimeSolution> solve(int samples, double t0, DoubleMatrix2D X0, double tf, boolean printProgress) throws Exception {
		return solve(samples, t0, X0, tf, 1.0, printProgress);
	}

	public List<FiniteTimeSolution> solve(int samples, final double t0, DoubleMatrix2D X0, final double tf, final double recordStep) throws Exception {
		return solve(samples, t0, X0, tf, recordStep, false);
	}

	public List<FiniteTimeSolution> solve(int samples, final double t0, DoubleMatrix2D X0, final double tf, final double recordStep, boolean printProgress) throws Exception {
		assert(samples > 0);
		ensurePoolStarted();
		ExecutorCompletionService<FiniteTimeSolution> ecs = new ExecutorCompletionService<>(executor);
		List<FiniteTimeSolution> solutions = new ArrayList<>(samples);
		for (int i=0; i < samples; i++) {
			final DoubleMatrix1D x0 = X0.viewColumn(i);
			final SSASolver solver = solverFactory.createSolver();
			Callable<FiniteTimeSolution> task = new Callable<FiniteTimeSolution>() {

				@Override
				public FiniteTimeSolution call() throws Exception {
					FiniteTimeSolution solution = solver.solve(t0, x0, tf, recordStep);
					return solution;
				}

			};
			ecs.submit(task);
		}
		for (int i=0; i < samples; i++) {
			Future<FiniteTimeSolution> future = ecs.take();
			FiniteTimeSolution solution = future.get();
			solutions.add(solution);
			if (printProgress) {
				System.out.print("\r" + (i+1) + " out of " + samples + " samples finished");
				System.out.flush();
			}
		}
		if (printProgress) {
			System.out.println();
		}
		return solutions;
	}

}
