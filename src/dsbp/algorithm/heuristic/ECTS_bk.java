package dsbp.algorithm.heuristic;

import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import dsbp.model.*;
import dsbp.util.*;

/**
 * Basic fixed size implementation for Tabu list
 * @author Hoang Tran
 *
 */


/**
 * This class is a Tabu Search implementation.
 *
 * @author Hoang Tran
 */
//public class ECTS_bk extends Heuristic {
//
//    /**
//     * Tabu parameters.
//     */
//	private int tabuListSize;
//	private TabuList tabuList;
//	/**
//     * Instantiates a new SA.
//     *
//     * @param problem problem reference
//     * @param random  random number generator.
//     * @param alpha   cooling rate for the simulated annealing
//     * @param t0      initial temperature, T0
//     * @param saMax   number of iterations before update the temperature
//     */
//    public ECTS_bk(Problem problem, Random random, int tabuListSize) {
//        super(problem, random, "TS");
//
//        // initializing simulated annealing parameters
//        this.tabuListSize = tabuListSize;
//        this.tabuList = new TabuList(this.tabuListSize);
//    }
//	
//    /**
//	 * Find the non-tabu {@link Solution} with the lowest value.<br>
//	 * This method doesn't use any Aspiration Criteria.
//	 */
//	
//    public Solution findBestNeighbor(List<Solution> neighborsSolutions, final List<Solution> solutionsInTabu) {
//		
//		//sort the neighbors
//		Collections.sort(neighborsSolutions, new Comparator<Solution>() {
//			@Override
//			public int compare(Solution a, Solution b) {
//				return ((Double)a.getTotalCost()).compareTo((Double)b.getTotalCost());
//			}
//		});
//
//    	//remove any neighbor that is in tabu list
//		//CollectionUtils.filterInverse(neighborsSolutions, new Predicate<Solution>() {
//		//	@Override
//		//	public boolean evaluate(Solution neighbor) {
//		//		if (solutionsInTabu.contains(neighbor)) {
//		//			System.out.println(count++);
//		//		}
//		//		return solutionsInTabu.contains(neighbor);
//		//	}
//		//});
//				
//		//get the neighbor with lowest value
//		int i = 0;
//		while (solutionsInTabu.contains(neighborsSolutions.get(i))) i++;
//		return neighborsSolutions.get(i);
//	}
//    /**
//     * Executes the Simulated Annealing.
//     *
//     * @param initialSolution the initial (input) solution.
//     * @param timeLimitMillis the time limit (in milliseconds).
//     * @param maxIters        the maximum number of iterations without improvements to execute.
//     * @param output          output PrintStream for logging purposes.
//     * @return the best solution encountered by the SA.
//     */
//    public Solution run(Solution initialSolution, long timeLimitMillis, long maxIters, PrintStream output) {
//        /*
//    	long finalTimeMillis = System.currentTimeMillis() + timeLimitMillis;
//
//        bestSolution = initialSolution;
//        Solution currentSolution = initialSolution.clone();
//        
//        int currentIteration = 0;
//        while (System.currentTimeMillis() < finalTimeMillis) {
//			//List<Solution> candidateNeighbors = currentSolution.getNeighbors();
//			List<Solution> candidateNeighbors = new LinkedList<Solution>();
//			List<Solution> solutionsInTabu = IteratorUtils.toList(tabuList.iterator());
//			
//			for (int i = 0; i < 500; i++) {
//				Solution neighbor = currentSolution.getRandomNeighbor(random);
//				candidateNeighbors.add(neighbor);
//			}
//			
//			Solution bestNeighborFound = findBestNeighbor(candidateNeighbors, solutionsInTabu);
//			if (bestNeighborFound.getTotalCost() < bestSolution.getTotalCost()) {
//				bestSolution = bestNeighborFound;
//				Util.safePrintStatus(output, nIters, bestSolution, currentSolution, "*");
//			}
//			tabuList.add(currentSolution);
//			currentSolution = bestNeighborFound;
//			tabuList.updateSize(currentIteration, bestSolution);
//            nIters++;
//        }
//		*/
//        return bestSolution;
//    }
//
//    /**
//     * Returns the string representation of this heuristic.
//     *
//     * @return the string representation of this heuristic (with parameters values).
//     */
//    public String toString() {
//        return String.format("Simulated Annealing (size=%d)", this.tabuListSize);
//    }
//}
