package dsbp.algorithm.heuristic;

import java.io.*;
import java.util.*;

import dsbp.model.*;

/**
 * This class is a Late Acceptance Hill Climbing implementation.
 *
 * @author Hoang Tran
 */
public class LAHC extends Heuristic {

	/**
	 * LAHC list with solution costs.
	 */
	private double[] list;


	/**
	 * Instantiates a new LAHC.
	 *
	 * @param problem  problem.
	 * @param random   random number generator.
	 * @param listSize LAHC list size.
	 */
	public LAHC(Problem problem, Random random, int listSize) {
		super(problem, random, "LAHC");

		// initializing the late acceptance list
		list = new double[listSize];
	}

	/**
	 * Executes the LAHC.
	 *
	 * @param initialSolution the initial (input) solution.
	 * @param timeLimitMillis the time limit (in milliseconds).
	 * @param maxIters        the maximum number of iterations without improvements to execute.
	 * @param output
	 * @return the best solution encountered by the LAHC.
	 */
	public Solution run(Solution initialSolution, long timeLimitMillis, long maxIters, PrintStream output) {
		long finalTimeMillis = System.currentTimeMillis() + timeLimitMillis;

		bestSolution = initialSolution;
		Solution solution = initialSolution.clone();

		// initializing LAHC list
		for (int i = 0; i < list.length; i++)
			list[i] = initialSolution.getTotalCost();

		int nItersWithoutImprovement = 0;
		int positionList = -1;

		//while (System.currentTimeMillis() < finalTimeMillis) {
		while (System.currentTimeMillis() < finalTimeMillis && nItersWithoutImprovement++ < maxIters) {
			positionList = (positionList + 1) % list.length;
			
			//Move move = selectMove(solution);
			//double delta = move.doMove(solution);
			Solution neighbor = solution.getRandomNeighbor(random);
			double delta = neighbor.getTotalCost() - solution.getTotalCost();
			// if solution is improved...
			if (delta < 0) {
				//acceptMove(move);
				solution = neighbor;
				nItersWithoutImprovement = 0;
				
				if (solution.getTotalCost() < bestSolution.getTotalCost()) {
					bestSolution = solution.clone();
					//Util.safePrintStatus(output, nIters, bestSolution, solution, "*");
				}
			}

			// if solution is not improved, but is accepted...
			else if (delta < 0.001 || neighbor.getTotalCost() <= list[positionList]) {
				//acceptMove(move);
				solution = neighbor;
			}

			// if solution is rejected..
			else {
				//rejectMove(move);
			}

			list[positionList] = solution.getTotalCost();
			nIters++;
		}

		//if (System.currentTimeMillis() < finalTimeMillis) {
		//	nItersWithoutImprovement = 0;
		//	for (int i = 0; i < list.length; i++)
		//		list[i] = initialSolution.getTotalCost();
		//	if (USE_LEARNING) learningAutomata.initProbabilities(getMoves());
		//	Util.safePrintText(output, "Resetting LAHC list", "");
		//}
		//}

		return bestSolution;
	}

	/**
	 * Returns the string representation of the heuristic.
	 *
	 * @return the string representation of the heuristic.
	 */
	public String toString() {
		return String.format("LAHC (listSize=%d)", list.length);
	}
}
