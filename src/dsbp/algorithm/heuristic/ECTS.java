package dsbp.algorithm.heuristic;

import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import dsbp.algorithm.neighborhood.ECTSMove;
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
public class ECTS extends Heuristic {

    /**
     * Tabu parameters.
     */
	private int tabuListSize;
	private HashSet<Integer> tabuSet;
	private LinkedList<Integer> tabuList;
	/**
     * Instantiates a new SA.
     *
     * @param problem problem reference
     * @param random  random number generator.
     */
    public ECTS(Problem problem, Random random) {
        super(problem, random, "ECTS");
        
        tabuSet = new HashSet<Integer>();
        tabuList = new LinkedList<Integer>();
    }
	
    /**
     * Executes the Tabu.
     *
     * @param initialSolution the initial (input) solution.
     * @param timeLimitMillis the time limit (in milliseconds).
     * @param maxIters        the maximum number of iterations without improvements to execute.
     * @param output          output PrintStream for logging purposes.
     * @return the best solution encountered by the TS.
     */
    public Solution run(Solution initialSolution, long timeLimitMillis, long maxIters, PrintStream output) {
        long finalTimeMillis = System.currentTimeMillis() + timeLimitMillis;
        ECTSMove move = new ECTSMove(problem, random, 1);
        
        bestSolution = initialSolution;
        Solution solution = initialSolution.clone();
        maxIters = 100 * solution.size; 
        tabuListSize = 10 * solution.size;		
        int nItersWithoutImprovement = 0;
        while (System.currentTimeMillis() < finalTimeMillis && nItersWithoutImprovement++ < maxIters) {
            move.doMove(solution, tabuList, tabuSet);
            acceptMove(move);
            if (solution.getTotalCost() < bestSolution.getTotalCost()) {
            	nItersWithoutImprovement = 0;
            	bestSolution = solution.clone();
            	Util.safePrintStatus(output, nIters, bestSolution, solution, "*");      
            }
            // updates tabu list
            tabuList.add(solution.code);
            tabuSet.add(solution.code);
            if (tabuList.size() > tabuListSize) {
            	tabuSet.remove(tabuList.get(0));
            	tabuList.remove(0);
            }
            
            nIters++;
        }
        /*
        int currentIteration = 0;
        while (System.currentTimeMillis() < finalTimeMillis) {
        	//List<Solution> candidateNeighbors = currentSolution.getNeighbors();
			List<Solution> candidateNeighbors = new LinkedList<Solution>();
			List<Solution> solutionsInTabu = IteratorUtils.toList(tabuList.iterator());
			
			for (int i = 0; i < 500; i++) {
				Solution neighbor = currentSolution.getRandomNeighbor(random);
				candidateNeighbors.add(neighbor);
			}
			
			Solution bestNeighborFound = findBestNeighbor(candidateNeighbors, solutionsInTabu);
			if (bestNeighborFound.getTotalCost() < bestSolution.getTotalCost()) {
				bestSolution = bestNeighborFound;
				Util.safePrintStatus(output, nIters, bestSolution, currentSolution, "*");
			}
			tabuList.add(currentSolution);
			currentSolution = bestNeighborFound;
			tabuList.updateSize(currentIteration, bestSolution);
            nIters++;
        }
		*/
        	
        return bestSolution;
    }

    /**
     * Returns the string representation of this heuristic.
     *
     * @return the string representation of this heuristic (with parameters values).
     */
    public String toString() {
        return String.format("Simulated Annealing (size=%d)", this.tabuListSize);
    }
}
