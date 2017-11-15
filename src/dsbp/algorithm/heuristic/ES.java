package dsbp.algorithm.heuristic;

import java.io.PrintStream;
import java.util.Random;

import dsbp.model.*;

/**
 * This class is a Evolutionary Strategy implementation.
 *
 * @author Hoang Tran
 */
public class ES extends Heuristic {

	public ES(Problem problem, Random random, String name) {
		super(problem, random, name);
		// TODO Auto-generated constructor stub
	}
//
//    /**
//     * ES parameters.
//     */
//	int[][][] solutionMapping; 
//    //CMAESOptimizer optimizer = new CMAESOptimizer(30, 0, true, 10, 
//    //		0, new MersenneTwister(), false, null);
//    
//	MultivariateFunction fitnessFunction = new MultivariateFunction() {
//        public double value(double[] dispatchSolution) {
//        	double drivingCost = 0;
//        	double balancingCost = 0;
//            
//        	// compute driving cost
//            for (int k = 0; k < problem.nTimeHorizons; k++) { 
//    			for (int i = 0; i < problem.nRegions; i++) {
//    				for (int j = 0; j < problem.nRegions; j++) {
//    					if (i < j && dispatchSolution[solutionMapping[k][i][j]] > 0) {
//    						drivingCost += dispatchSolution[solutionMapping[k][i][j]] * problem.distances[i][j];
//    					} else if (j < i && dispatchSolution[solutionMapping[k][j][i]] < 0) {
//    						drivingCost += -dispatchSolution[solutionMapping[k][j][i]] * problem.distances[i][j];
//    					}			
//    				}
//    			}
//    		}
//            
//            // compute balancing cost
//            int[] supply = new int[problem.nRegions];
//            for (int i = 0; i < problem.nRegions; i++) {
//            	supply[i] = problem.initialSupply[i];
//            }
//            for (int k = 0; k < problem.nTimeHorizons; k++) {
//            	// update supply after dispatch
//            	for (int i = 0; i < problem.nRegions; i++) {
//    				for (int j = 0; j < problem.nRegions; j++) {
//    					if (i < j && dispatchSolution[solutionMapping[k][i][j]] > 0) {
//    						supply[i] -= dispatchSolution[solutionMapping[k][i][j]];
//    						supply[j] += dispatchSolution[solutionMapping[k][i][j]];
//    					} else if (j < i && dispatchSolution[solutionMapping[k][j][i]] < 0) {
//    						supply[i] -= -dispatchSolution[solutionMapping[k][j][i]];
//    						supply[j] += -dispatchSolution[solutionMapping[k][j][i]];
//    					}
//    				}
//    			}
//            	// compute cost
//    			for (int i = 0; i < problem.nRegions; i++) {
//    				if (supply[i] > 0) 
//    					balancingCost += Math.abs(problem.predictedDemand[k][i] * 1.0 / supply[i] - problem.getTotalDemand(k) * 1.0 / problem.getTotalSupply());
//    				else
//    					balancingCost += Math.abs(problem.predictedDemand[k][i] * 1.0 / 1 - problem.getTotalDemand(k) * 1.0 / problem.getTotalSupply());
//    			}
//    			// update supply after mobility
//    			double[] temp = new double[problem.nRegions];
//    			for (int ii = 0; ii < problem.nRegions; ii++) {
//    				for (int jj = 0; jj < problem.nRegions; jj++) {
//    					temp[ii] += supply[jj] * problem.mobilityPattern[jj][ii];
//    				}
//    			}
//    			for (int ii = 0; ii < problem.nRegions; ii++) {
//    				supply[ii] = (int) Math.round(temp[ii]);
//    			}
//    		}
//    		
//            double totalCost = drivingCost + problem.weightFactor * balancingCost; 
//            return totalCost;
//        }
//    };	
//	
//    /**
//     * Instantiates a new ES.
//     *
//     * @param problem problem reference
//     * @param random  random number generator.
//     */
//    public ES(Problem problem, Random random) {
//        super(problem, random, "ES");
//
//        // initializing evolutionary strategy parameters
//    	this.solutionMapping = new int[problem.nTimeHorizons][problem.nRegions][problem.nRegions];
//    	for (int k = 0; k < problem.nTimeHorizons; k++) { 
//    		for (int i = 0; i < problem.nRegions - 1; i++) {
//    			for (int j = i + 1; j < problem.nRegions; j++) {
//    				this.solutionMapping[k][i][j] = k * problem.nRegions * (problem.nRegions - 1) / 2 + i * problem.nRegions - i * (i + 1) /2 + j - i - 1;
//    			}
//    		}
//    	}
//
//    }
//
//	
//     /**
//     * Executes the Evolutionary Strategy.
//     *
//     * @param initialSolution the initial (input) solution.
//     * @param timeLimitMillis the time limit (in milliseconds).
//     * @param maxIters        the maximum number of iterations without improvements to execute.
//     * @param output          output PrintStream for logging purposes.
//     * @return the best solution encountered by the SA.
//     */
//    public Solution run(Solution initialSolution, long timeLimitMillis, long maxIters, PrintStream output) {
//        //long finalTimeMillis = System.currentTimeMillis() + timeLimitMillis;
//
//        bestSolution = initialSolution;
//	    double[] start = new double[initialSolution.size];
//	    double[] lower = new double[initialSolution.size];
//	    double[] upper = new double[initialSolution.size];
//	    double[] sigma = new double[initialSolution.size];
//	    for (int i = 0; i < initialSolution.size; i++) {
//	    	start[i] = 0;
//	    	lower[i] = -10;
//	    	upper[i] = 10;
//	    	sigma[i] = 0.1;
//	    }
//	    double[] result = optimizer.optimize(new MaxEval(10),
//	                                               new ObjectiveFunction(fitnessFunction),
//	                                               GoalType.MINIMIZE,
//	                                               new CMAESOptimizer.PopulationSize(5),
//	                                               new CMAESOptimizer.Sigma(sigma),
//	                                               new InitialGuess(start),
//	                                               new SimpleBounds(lower, upper)).getPoint();
//	    for (int i = 0; i < initialSolution.size; i++) {
//	    	bestSolution.dispatch[i] = (int) Math.round(result[i]);
//	    }
//	    System.out.println(Boolean.toString(bestSolution.validate()));
//        return bestSolution;
//    }
//
//    /**
//     * Returns the string representation of this heuristic.
//     *
//     * @return the string representation of this heuristic (with parameters values).
//     */
//    public String toString() {
//        return String.format("Evolutionary Strategy (size=)");
//    }

	@Override
	public Solution run(Solution solution, long timeLimitMillis, long maxIters, PrintStream output) {
		// TODO Auto-generated method stub
		return null;
	}

}