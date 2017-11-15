package dsbp.algorithm.heuristic;

import java.io.*;
import java.util.*;

import dsbp.model.*;
import dsbp.util.*;

/**
 * This class is a Step Counting Hill Climbing implementation.
 *
 * @author Tulio Toffolo
 */
public class SCHC extends Heuristic {

    /**
     * SCHC parameters.
     */
    private int stepSize;
    private double costBound;

    /**
     * Instantiates a new SCHC.
     *
     * @param problem  problem
     * @param random   random number generator
     * @param stepSize the SCHC step size value
     */
    public SCHC(Problem problem, Random random, int stepSize ) {
        super(problem, random, "SCHC");

        // initializing the late acceptance list
        this.stepSize = stepSize;
        this.costBound = Double.MAX_VALUE;
    }

    /**
     * Executes the SCHC.
     *
     * @param initialSolution the initial (input) solution.
     * @param timeLimitMillis the time limit (in milliseconds).
     * @param maxIters        the maximum number of iterations without improvements to execute.
     * @param output
     * @return the best solution encountered by the SCHC.
     */
    public Solution run(Solution initialSolution, long timeLimitMillis, long maxIters, PrintStream output) {
        long finalTimeMillis = System.currentTimeMillis() + timeLimitMillis;

        bestSolution = initialSolution;
        Solution solution = initialSolution.clone();

        costBound = bestSolution.getTotalCost();

        int nItersWithoutImprovement = 0;
        int stepCounter = 0;

        //while (System.currentTimeMillis() < finalTimeMillis) {
        while (System.currentTimeMillis() < finalTimeMillis && nItersWithoutImprovement++ < maxIters) {
            stepCounter++;

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
            else if (delta < 0.001 || neighbor.getTotalCost() <= costBound) {
                //acceptMove(move);
            	solution = neighbor;
            }

            // if solution is rejected..
            else {
                //rejectMove(move);
            }

            // if necessary, updates costBound
            if (stepCounter >= stepSize) {
                costBound = solution.getTotalCost();
                stepCounter = 0;
            }

            nIters++;
        }

        //if (System.currentTimeMillis() < finalTimeMillis) {
        //    nItersWithoutImprovement = 0;
        //    costBound = initialSolution.getTotalCost();
        //    stepCounter = 0;
        //    Util.safePrintText(output, "Restarting SCHC cost bound", "");
        //}
        //}

        return bestSolution;
    }


    /**
     * Returns the string representation of the heuristic.
     *
     * @return the string representation of the heuristic (with parameters).
     */
    public String toString() {
        return String.format("SCHC (stepSize=%s)", Util.longToString(stepSize));
    }
}
