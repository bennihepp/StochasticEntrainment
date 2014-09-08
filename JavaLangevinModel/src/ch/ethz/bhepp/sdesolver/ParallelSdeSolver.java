package ch.ethz.bhepp.sdesolver;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class ParallelSdeSolver {

	private static final int DEFAULT_NUM_OF_THREADS = 4;
	private int numOfThreads;
	private ExecutorService executor;
	private SdeStepperFactory stepperFactory;

	public ParallelSdeSolver(SdeStepperFactory stepperFactory) {
		this(stepperFactory, 0);
	}

	public ParallelSdeSolver(SdeStepperFactory stepperFactory, int numOfThreads) {
		this.stepperFactory = stepperFactory;
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

	public List<SdeSolution> solve(int samples, double t0, DoubleMatrix2D X0, double tf) throws Exception {
		return solve(samples, t0, X0, tf, 1.0);
	}

	public List<SdeSolution> solve(int samples, final double t0, DoubleMatrix2D X0, final double tf, final double recordStep) throws Exception {
		assert(samples > 0);
		ensurePoolStarted();
		ExecutorCompletionService<SdeSolution> ecs = new ExecutorCompletionService<>(executor);
		List<SdeSolution> solutions = new ArrayList<>(samples);
		for (int i=0; i < samples; i++) {
			SdeStepper stepper = stepperFactory.createStepper();
			final DoubleMatrix1D x0 = X0.viewColumn(i);
			final SdeSolver solver = new SdeSolver(stepper);
			Callable<SdeSolution> task = new Callable<SdeSolution>() {

				@Override
				public SdeSolution call() throws Exception {
					SdeSolution solution = solver.solve(t0, x0, tf, recordStep);
					return solution;
				}

			};
			ecs.submit(task);
		}
		for (int i=0; i < samples; i++) {
			Future<SdeSolution> future = ecs.take();
			SdeSolution solution = future.get();
			solutions.add(solution);
		}
		return solutions;
	}

}
