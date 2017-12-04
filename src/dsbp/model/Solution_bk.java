package dsbp.model;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import dsbp.util.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.lang.Math;

/**
 * This class represents a Solution of the Unrelated Parallel Machine Scheduling
 * Problem..
 *
 * @author Hoang Tran
 */
public class Solution_bk {

	public final Problem problem;
	public final int[] dispatchSolution;
	public final int solutionSize;
	protected double drivingCost;
	protected double balancingCost;
	protected double totalCost;
	private int[][][] matrixToVectorMapping;
	private int[][] vectorToMatrixMapping;
	private int[][] dispatchSupply;
	private double[][] demandSupplyRatio;
	private int hashCode;

	/**
	 * Instantiates a new Solution.
	 *
	 * @param problem problem considered.
	 */
	public Solution_bk(Problem problem) {
		this.problem = problem;
		this.solutionSize = problem.nTimeHorizons * problem.nRegions * (problem.nRegions - 1) / 2;
		this.dispatchSolution = new int[solutionSize];
		this.drivingCost = 0;
		this.balancingCost = 0;
		this.totalCost = 0;
		this.matrixToVectorMapping = new int[problem.nTimeHorizons][problem.nRegions][problem.nRegions];
		this.vectorToMatrixMapping = new int[solutionSize][3];
		for (int k = 0; k < problem.nTimeHorizons; k++) { 
			for (int i = 0; i < problem.nRegions - 1; i++) {
				for (int j = i + 1; j < problem.nRegions; j++) {
					this.matrixToVectorMapping[k][i][j] = k * problem.nRegions * (problem.nRegions - 1) / 2 + i * problem.nRegions - i * (i + 1) /2 + j - i - 1;
					this.vectorToMatrixMapping[this.matrixToVectorMapping[k][i][j]][0] = k;
					this.vectorToMatrixMapping[this.matrixToVectorMapping[k][i][j]][1] = i;
					this.vectorToMatrixMapping[this.matrixToVectorMapping[k][i][j]][2] = j;
				}
			}
		}
		this.dispatchSupply = new int[problem.nTimeHorizons][problem.nRegions];
		this.demandSupplyRatio = new double[problem.nTimeHorizons][problem.nRegions];
		this.hashCode = Arrays.hashCode(this.dispatchSolution);
		updateTotalCost();
	}

	/**
	 * Private constructor used for cloning.
	 *
	 * @param solution solution to copy from.
	 */
	private Solution_bk(Solution_bk solution) {
		this.problem = solution.problem;
		this.solutionSize = solution.solutionSize;
		this.dispatchSolution = new int[solutionSize];
		for (int i = 0; i < solutionSize; i++) {
			this.dispatchSolution[i] = solution.dispatchSolution[i];
		}
		this.drivingCost = solution.drivingCost;
		this.balancingCost = solution.balancingCost;
		this.totalCost = solution.totalCost;
		this.matrixToVectorMapping = solution.matrixToVectorMapping;
		this.vectorToMatrixMapping = solution.vectorToMatrixMapping;
		this.dispatchSupply = new int[problem.nTimeHorizons][problem.nRegions];
		for (int k = 0; k < problem.nTimeHorizons; k++) {
			for (int i = 0; i < problem.nRegions; i++) {
				this.dispatchSupply[k][i] = solution.dispatchSupply[k][i];
			}
		}
		this.demandSupplyRatio = new double[problem.nTimeHorizons][problem.nRegions];
		for (int k = 0; k < problem.nTimeHorizons; k++) {
			for (int i = 0; i < problem.nRegions; i++) {
				this.demandSupplyRatio[k][i] = solution.demandSupplyRatio[k][i];
			}
		}
		this.hashCode = solution.hashCode;
	}

	/**
	 * Creates and returns a copy of this solution.
	 */
	public Solution_bk clone() {
		return new Solution_bk(this);
	}

	/**
	 * Creates and returns a random neighbor of this solution.
	 */
	public Solution_bk getRandomNeighbor(Random random) {
		int position;
		int alpha;
		int num = random.nextInt(3) + 1;
		do {
			position = random.nextInt(solutionSize);
			alpha = random.nextInt(2)==0 ? num : -num;
		} while (!validate(position, alpha));
		Solution_bk neighbor = new Solution_bk(this);
		neighbor.dispatchSolution[position] += alpha;
		neighbor.updateTotalCost();
		return neighbor;
	}

	/**
	 * Do a random move of this solution.
	 */
	public void doRandomMove(Random random) {
		int position;
		int alpha;
		boolean valid;
		do {
			valid = false;
			position = random.nextInt(solutionSize);
			int k = vectorToMatrixMapping[position][0];
			int i = vectorToMatrixMapping[position][1];
			int j = vectorToMatrixMapping[position][2];
			int num = random.nextInt(3) + 1;
			alpha = demandSupplyRatio[k][i] < demandSupplyRatio[k][j] ? num : -num;
			if (alpha > 0 && dispatchSupply[k][i] - alpha >= 0) valid = true;
			else if (alpha < 0 && dispatchSupply[k][j] + alpha >= 0) valid = true;
		} while (!valid);
		this.dispatchSolution[position] += alpha;
		this.updateTotalCost();
	}

	/**
	 * Gets the solution driving cost. 
	 *
	 * @return the solution driving cost.
	 */
	public double getDrivingCost() {
		return drivingCost;
	}

	/**
	 * Gets the solution balancing cost. 
	 *
	 * @return the solution balancing cost.
	 */
	public double getBalancingCost() {
		return balancingCost;
	}

	/**
	 * Gets the solution total cost. 
	 *
	 * @return the solution total cost.
	 */
	public double getTotalCost() {
		return totalCost;
	}

	/**
	 * Updates (and returns) the driving cost of the solution.
	 *
	 * @return the updated solution cost.
	 */
	private double updateDrivingCost() {
		// compute driving cost
		drivingCost = 0;
		for (int k = 0; k < problem.nTimeHorizons; k++) { 
			for (int i = 0; i < problem.nRegions; i++) {
				for (int j = 0; j < problem.nRegions; j++) {
					if (i < j && dispatchSolution[matrixToVectorMapping[k][i][j]] > 0) {
						drivingCost += dispatchSolution[matrixToVectorMapping[k][i][j]] * problem.distances[i][j];
					} else if (j < i && dispatchSolution[matrixToVectorMapping[k][j][i]] < 0) {
						drivingCost += -dispatchSolution[matrixToVectorMapping[k][j][i]] * problem.distances[i][j];
					}			
				}
			}
		}
		return drivingCost;
	}

	/**
	 * Updates (and returns) the balancing cost of the solution.
	 *
	 * @return the updated solution cost.
	 */
	private double updateBalancingCost() {
		// compute balancing cost
		balancingCost = 0;
		int[] supply = new int[problem.nRegions];
		for (int i = 0; i < problem.nRegions; i++) {
			supply[i] = problem.initialSupply[i];
		}
		for (int k = 0; k < problem.nTimeHorizons; k++) {
			// update supply after dispatch
			for (int i = 0; i < problem.nRegions; i++) {
				for (int j = 0; j < problem.nRegions; j++) {
					if (i < j && dispatchSolution[matrixToVectorMapping[k][i][j]] > 0) {
						supply[i] -= dispatchSolution[matrixToVectorMapping[k][i][j]];
						supply[j] += dispatchSolution[matrixToVectorMapping[k][i][j]];
					} else if (j < i && dispatchSolution[matrixToVectorMapping[k][j][i]] < 0) {
						supply[i] -= -dispatchSolution[matrixToVectorMapping[k][j][i]];
						supply[j] += -dispatchSolution[matrixToVectorMapping[k][j][i]];
					}
				}
			}
			// compute cost
			for (int i = 0; i < problem.nRegions; i++) {
				dispatchSupply[k][i] = supply[i];
				if (supply[i] > 0)
					demandSupplyRatio[k][i] = problem.predictedDemand[k][i] * 1.0 / supply[i];
				else
					demandSupplyRatio[k][i] = problem.predictedDemand[k][i] * 1.0 / 1;
				balancingCost += Math.abs(demandSupplyRatio[k][i] - problem.getGlobalDemandSupplyRatio(k));
			}
			// update supply after mobility
			double[] temp = new double[problem.nRegions];
			for (int ii = 0; ii < problem.nRegions; ii++) {
				for (int jj = 0; jj < problem.nRegions; jj++) {
					temp[ii] += supply[jj] * problem.mobilityPattern[jj][ii];
				}
			}
			for (int ii = 0; ii < problem.nRegions; ii++) {
				supply[ii] = (int) Math.round(temp[ii]);
			}
		}
		return balancingCost;
	}

	/**
	 * Updates (and returns) the total cost of the solution.
	 *
	 * @return the updated solution cost.
	 */
	public double updateTotalCost() {
		updateDrivingCost();
		updateBalancingCost();
		totalCost = drivingCost + problem.weightFactor * balancingCost;
		hashCode = Arrays.hashCode(dispatchSolution);
		return totalCost;
	}

	public void updateTotalCost(double totalCost) {
		this.totalCost = totalCost;
		
	}
	/**
	 * Validates the solution.
	 *
	 * @param output output stream (example: System.out) to print eventual error
	 *               messages.
	 * @return true if the solution and its costs are valid and false otherwise.
	 */
	public boolean validate() {
		boolean valid = true;

		// initial supply
		int[] supply = new int[problem.nRegions];
		for (int i = 0; i < problem.nRegions; i++) {
			supply[i] = problem.initialSupply[i];
		}

		for (int k = 0; k < problem.nTimeHorizons; k++) { 
			// update supply after dispatch
			for (int i = 0; i < problem.nRegions; i++) {
				for (int j = 0; j < problem.nRegions; j++) {
					if (i < j && dispatchSolution[matrixToVectorMapping[k][i][j]] > 0) {
						supply[i] -= dispatchSolution[matrixToVectorMapping[k][i][j]];
						supply[j] += dispatchSolution[matrixToVectorMapping[k][i][j]];
					} else if (j < i && dispatchSolution[matrixToVectorMapping[k][j][i]] < 0) {
						supply[i] -= -dispatchSolution[matrixToVectorMapping[k][j][i]];
						supply[j] += -dispatchSolution[matrixToVectorMapping[k][j][i]];
					}
				}
			}
			// update supply after dispatch
			for (int i = 0; i < problem.nRegions; i++) {
				if (supply[i] < 0) {
					valid = false;
					//Util.safePrintf(output, "Supply %d at time %d is nagative\n", i, k);
				}
			}

			// update supply after mobility
			double[] temp = new double[problem.nRegions];
			for (int ii = 0; ii < problem.nRegions; ii++) {
				for (int jj = 0; jj < problem.nRegions; jj++) {
					temp[ii] += supply[jj] * problem.mobilityPattern[jj][ii];
				}
			}
			for (int ii = 0; ii < problem.nRegions; ii++) {
				supply[ii] = (int) Math.round(temp[ii]);
			}
		}

		return valid;
	}

	/**
	 * Validates the solution.
	 *
	 * @param output output stream (example: System.out) to print eventual error
	 *               messages.
	 * @return true if the solution and its costs are valid and false otherwise.
	 */
	public boolean validate(int position, int alpha) {
		boolean valid = false;
		int k = vectorToMatrixMapping[position][0];
		int i = vectorToMatrixMapping[position][1];
		int j = vectorToMatrixMapping[position][2];
		if (alpha > 0 && dispatchSupply[k][i] - alpha >= 0) valid = true;
		else if (alpha < 0 && dispatchSupply[k][j] + alpha >= 0) valid = true;
		return valid;
	}

	public List<Solution_bk> getNeighbors() {
		int alpha = 0;
		List<Solution_bk> neighbors = new LinkedList<Solution_bk>();
		for (int i = 0; i < this.solutionSize; i++) {
			Solution_bk neighbor1 = this.clone();
			Solution_bk neighbor2 = this.clone();
			neighbor1.dispatchSolution[i] += alpha;
			if (neighbor1.validate()) {
				neighbor1.updateTotalCost();
				neighbors.add(neighbor1);
			}
			neighbor2.dispatchSolution[i] -= alpha;
			if (neighbor2.validate()) {
				neighbor2.updateTotalCost();
				neighbors.add(neighbor2);
			}
		}
		return neighbors;
	}

	@Override
	public boolean equals(Object o) {
		return ((Integer) this.hashCode()).equals(((Solution_bk)o).hashCode());
	}

	@Override
	public int hashCode() {
		return hashCode;
	}

	/**
	 * Reads the solution from a file.
	 *
	 * @param filePath the output file path.
	 * @throws IOException in case any IO error occurs.
	 */
	public void read(String filePath) throws IOException {
		BufferedReader reader = Files.newBufferedReader(Paths.get(filePath));
		SimpleTokenizer token;
		for (int k = 0; k < problem.nTimeHorizons; k++) { 
			for (int i = 0; i < problem.nRegions; i++) {
				token = new SimpleTokenizer(reader.readLine());
				for (int j = 0; j < problem.nRegions; j++) {
					int h = token.nextInteger();
					if (i < j && h > 0) {
						dispatchSolution[matrixToVectorMapping[k][i][j]] = h;
					} else if (j < i && h > 0) {
						dispatchSolution[matrixToVectorMapping[k][j][i]] = -h;
					}
				}
			}
			reader.readLine();
		}
		updateTotalCost();
		reader.close();
	}

	/**
	 * Writes the solution to a file.
	 *
	 * @param filePath the output file path.
	 * @throws IOException in case any IO error occurs.
	 */
	public void write(String filePath) throws IOException {
		PrintWriter writer = new PrintWriter(Files.newBufferedWriter(Paths.get(filePath)));

		// print dispatch solution
		writer.printf("Dispatch solution\n");
		for (int i = 0; i < problem.nTimeHorizons * problem.nRegions * (problem.nRegions - 1) / 2; i++) { 
			writer.printf("%d ", dispatchSolution[i]);
		}
		writer.printf("\n\n");

		// print supply
		int[] supply = new int[problem.nRegions];
		for (int i = 0; i < problem.nRegions; i++) {
			supply[i] = problem.initialSupply[i];
		}
		for (int k = 0; k < problem.nTimeHorizons; k++) { 
			writer.printf("k = %d\n", k+1);
			// print dispatch solution
			/*
        	writer.printf("Dispatch solution\n");
			for (int i = 0; i < problem.nRegions; i++) {
				for (int j = 0; j < problem.nRegions; j++) {
					if (j < i) {
						if (dispatchSolution[solutionMapping[k][j][i]] < 0) {
							writer.printf("%d ", -dispatchSolution[solutionMapping[k][j][i]]);
						} else {
							writer.printf("0 ");
						}
					} else if (j == i) {
						writer.printf("0 ");
					} else if (j > i) {
						if (dispatchSolution[solutionMapping[k][i][j]] > 0) {
							writer.printf("%d ", dispatchSolution[solutionMapping[k][i][j]]);
						} else {
							writer.printf("0 ");
						}
					}
				}
				writer.printf("\n");
			}
			 */
			// print demand
			writer.printf("Demand   ");
			for (int i = 0; i < problem.nRegions; i++) {
				writer.printf("%d ", problem.predictedDemand[k][i]);
			}
			// print supply after dispatch
			for (int i = 0; i < problem.nRegions; i++) {
				for (int j = 0; j < problem.nRegions; j++) {
					if (i < j && dispatchSolution[matrixToVectorMapping[k][i][j]] > 0) {
						supply[i] -= dispatchSolution[matrixToVectorMapping[k][i][j]];
						supply[j] += dispatchSolution[matrixToVectorMapping[k][i][j]];
					} else if (j < i && dispatchSolution[matrixToVectorMapping[k][j][i]] < 0) {
						supply[i] -= -dispatchSolution[matrixToVectorMapping[k][j][i]];
						supply[j] += -dispatchSolution[matrixToVectorMapping[k][j][i]];
					}
				}
			}
			writer.printf("\nSupply   ");
			for (int i = 0; i < problem.nRegions; i++) {
				writer.printf("%d ", supply[i]);
			}
			// print supply after mobility
			writer.printf("\nMobility ");
			double[] temp = new double[problem.nRegions];
			for (int ii = 0; ii < problem.nRegions; ii++) {
				for (int jj = 0; jj < problem.nRegions; jj++) {
					temp[ii] += supply[jj] * problem.mobilityPattern[jj][ii];
				}
			}
			for (int ii = 0; ii < problem.nRegions; ii++) {
				supply[ii] = (int) Math.round(temp[ii]);
				writer.printf("%d ", supply[ii]);
			}
			writer.printf("\n\n");
		}

		updateTotalCost();
		writer.printf("Driving cost=%f\n", drivingCost);
		writer.printf("Balancing cost=%f\n", balancingCost);
		writer.printf("Total cost=%f\n", totalCost);
		writer.close();
	}
}

// backup

///**
// * This class represents a Solution of the Unrelated Parallel Machine Scheduling
// * Problem..
// *
// * @author Hoang Tran
// */
//public class Solution {
//
//	public final Problem problem;
//    public final int[] dispatchSolution;
//    public final int solutionSize;
//    protected Double JDCost;
//    protected Double JECost;
//    protected Double dispatchCost;
//    private int[][][] solutionMapping;
//	
//    /**
//     * Instantiates a new Solution.
//     *
//     * @param problem problem considered.
//     */
//    public Solution(Problem problem) {
//        this.problem = problem;
//
//        solutionSize = problem.nTimeHorizons * problem.nRegions * (problem.nRegions - 1) / 2;
//        dispatchSolution = new int[solutionSize];
//		dispatchCost = 0.0;
//		solutionMapping = new int[problem.nTimeHorizons][problem.nRegions][problem.nRegions];
//		for (int k = 0; k < problem.nTimeHorizons; k++) { 
//			for (int i = 0; i < problem.nRegions - 1; i++) {
//				for (int j = i + 1; j < problem.nRegions; j++) {
//					solutionMapping[k][i][j] = k * problem.nRegions * (problem.nRegions - 1) / 2 + i * problem.nRegions - i * (i + 1) /2 + j - i - 1;
//				}
//			}
//		}
//    }
//
//    /**
//     * Private constructor used for cloning.
//     *
//     * @param solution solution to copy from.
//     */
//    private Solution(Solution solution) {
//        this.problem = solution.problem;
//
//        solutionSize = problem.nTimeHorizons * problem.nRegions * (problem.nRegions - 1) / 2;
//        dispatchSolution = new int[solutionSize];
//		for (int i = 0; i < solutionSize; i++) {
//            dispatchSolution[i] = solution.dispatchSolution[i];
//        }
//		solutionMapping = new int[problem.nTimeHorizons][problem.nRegions][problem.nRegions];
//		for (int k = 0; k < problem.nTimeHorizons; k++) { 
//			for (int i = 0; i < problem.nRegions - 1; i++) {
//				for (int j = i + 1; j < problem.nRegions; j++) {
//					solutionMapping[k][i][j] = k * problem.nRegions * (problem.nRegions - 1) / 2 + i * problem.nRegions - i * (i + 1) /2 + j - i - 1;
//				}
//			}
//		}
//		this.dispatchCost = solution.dispatchCost;
//    }
//
//    /**
//     * Creates and returns a copy of this solution.
//     */
//    public Solution clone() {
//        return new Solution(this);
//    }
//
//    public void readFromFile(String instancePath) throws IOException {
//        BufferedReader reader = Files.newBufferedReader(Paths.get(instancePath));
//
//        SimpleTokenizer token;
//        
//        for (int k = 0; k < problem.nTimeHorizons; k++) { 
//			for (int i = 0; i < problem.nRegions; i++) {
//				token = new SimpleTokenizer(reader.readLine());
//				for (int j = 0; j < problem.nRegions; j++) {
//					int h = token.nextInteger();
//					if (i < j && h > 0) {
//						dispatchSolution[solutionMapping[k][i][j]] = h;
//					} else if (j < i && h > 0) {
//						dispatchSolution[solutionMapping[k][j][i]] = -h;
//					}
//				}
//			}
//	        reader.readLine();
//		}
//        
//        updateCost();
//        
//        reader.close();
//    }
//    
//    /**
//     * Gets the solution makespan. Note the the makespan may be outdated if the
//     * solution was modified. To ensure that it is updated, call {@link
//     * #updateCost()}.
//     *
//     * @return the solution cost.
//     */
//    public Double getCost() {
//        return dispatchCost;
//    }
//    
//    public Double getJD() {
//        return JDCost;
//    }
//
//    public Double getJE() {
//        return JECost;
//    }
//    
//    /**
//     * Updates (and returns) the makespan of the solution.
//     *
//     * @return the updated solution cost.
//     */
//    public double updateCost() {
//    	double drivingCost = 0;
//    	double balancingCost = 0;
//        
//    	// compute driving cost
//        for (int k = 0; k < problem.nTimeHorizons; k++) { 
//			for (int i = 0; i < problem.nRegions; i++) {
//				for (int j = 0; j < problem.nRegions; j++) {
//					if (i < j && dispatchSolution[solutionMapping[k][i][j]] > 0) {
//						drivingCost += dispatchSolution[solutionMapping[k][i][j]] * problem.distances[i][j];
//					} else if (j < i && dispatchSolution[solutionMapping[k][j][i]] < 0) {
//						drivingCost += -dispatchSolution[solutionMapping[k][j][i]] * problem.distances[i][j];
//					}			
//				}
//			}
//		}
//        
//        // compute balancing cost
//        int[] supply = new int[problem.nRegions];
//        for (int i = 0; i < problem.nRegions; i++) {
//        	supply[i] = problem.initialSupply[i];
//        }
//        for (int k = 0; k < problem.nTimeHorizons; k++) {
//        	// update supply after dispatch
//        	for (int i = 0; i < problem.nRegions; i++) {
//				for (int j = 0; j < problem.nRegions; j++) {
//					if (i < j && dispatchSolution[solutionMapping[k][i][j]] > 0) {
//						supply[i] -= dispatchSolution[solutionMapping[k][i][j]];
//						supply[j] += dispatchSolution[solutionMapping[k][i][j]];
//					} else if (j < i && dispatchSolution[solutionMapping[k][j][i]] < 0) {
//						supply[i] -= -dispatchSolution[solutionMapping[k][j][i]];
//						supply[j] += -dispatchSolution[solutionMapping[k][j][i]];
//					}
//				}
//			}
//        	// compute cost
//			for (int i = 0; i < problem.nRegions; i++) {
//				if (supply[i] > 0) 
//					balancingCost += Math.abs(problem.predictedDemand[k][i] * 1.0 / supply[i] - problem.getTotalDemand(k) * 1.0 / problem.getTotalSupply());
//				else
//					balancingCost += Math.abs(problem.predictedDemand[k][i] * 1.0 / 1 - problem.getTotalDemand(k) * 1.0 / problem.getTotalSupply());
//			}
//			// update supply after mobility
//			double[] temp = new double[problem.nRegions];
//			for (int ii = 0; ii < problem.nRegions; ii++) {
//				for (int jj = 0; jj < problem.nRegions; jj++) {
//					temp[ii] += supply[jj] * problem.mobilityPattern[jj][ii];
//				}
//			}
//			for (int ii = 0; ii < problem.nRegions; ii++) {
//				supply[ii] = (int) Math.round(temp[ii]);
//			}
//		}
//		
//        //JDCost = Math.round(drivingCost);
//        //JECost = Math.round(balancingCost);
//        //dispatchCost = Math.round(drivingCost + problem.weightFactor * balancingCost); 
//        JDCost = drivingCost;
//        JECost = balancingCost;
//        dispatchCost = drivingCost + problem.weightFactor * balancingCost; 
//        return dispatchCost;
//    }
//
//    /**
//     * Validates the solution.
//     *
//     * @param output output stream (example: System.out) to print eventual error
//     *               messages.
//     * @return true if the solution and its costs are valid and false otherwise.
//     */
//    public boolean validate() {
//        boolean valid = true;
//        
//        // initial supply
//		int[] supply = new int[problem.nRegions];
//		for (int i = 0; i < problem.nRegions; i++) {
//			supply[i] = problem.initialSupply[i];
//		}
//		
//		for (int k = 0; k < problem.nTimeHorizons; k++) { 
//			// update supply after dispatch
//        	for (int i = 0; i < problem.nRegions; i++) {
//				for (int j = 0; j < problem.nRegions; j++) {
//					if (i < j && dispatchSolution[solutionMapping[k][i][j]] > 0) {
//						supply[i] -= dispatchSolution[solutionMapping[k][i][j]];
//						supply[j] += dispatchSolution[solutionMapping[k][i][j]];
//					} else if (j < i && dispatchSolution[solutionMapping[k][j][i]] < 0) {
//						supply[i] -= -dispatchSolution[solutionMapping[k][j][i]];
//						supply[j] += -dispatchSolution[solutionMapping[k][j][i]];
//					}
//				}
//			}
//			// update supply after dispatch
//			for (int i = 0; i < problem.nRegions; i++) {
//				if (supply[i] < 0) {
//					valid = false;
//	                //Util.safePrintf(output, "Supply %d at time %d is nagative\n", i, k);
//				}
//			}
//			
//			// update supply after mobility
//			double[] temp = new double[problem.nRegions];
//			for (int ii = 0; ii < problem.nRegions; ii++) {
//				for (int jj = 0; jj < problem.nRegions; jj++) {
//					temp[ii] += supply[jj] * problem.mobilityPattern[jj][ii];
//				}
//			}
//			for (int ii = 0; ii < problem.nRegions; ii++) {
//				supply[ii] = (int) Math.round(temp[ii]);
//			}
//		}
//
//        return valid;
//    }
//
//    public List<Solution> getNeighbors() {
//    	int alpha = 0;
//    	List<Solution> neighbors = new LinkedList<Solution>();
//    	for (int i = 0; i < this.solutionSize; i++) {
//    		Solution neighbor1 = this.clone();
//    		Solution neighbor2 = this.clone();
//    		neighbor1.dispatchSolution[i] += alpha;
//    		if (neighbor1.validate()) {
//    			neighbor1.updateCost();
//    			neighbors.add(neighbor1);
//    		}
//    		neighbor2.dispatchSolution[i] -= alpha;
//    		if (neighbor2.validate()) {
//    			neighbor2.updateCost();
//    			neighbors.add(neighbor2);
//    		}
//    	}
//    	return neighbors;
//    }
//    
//    @Override
//    public boolean equals(Object o) {
//    	return Arrays.equals(this.dispatchSolution, ((Solution)o).dispatchSolution);
//    }
//    
//    @Override
//    public int hashCode() {
//    	return this.dispatchSolution.hashCode();
//    }
//    
//    /**
//     * Writes the solution to a file.
//     *
//     * @param filePath the output file path.
//     * @throws IOException in case any IO error occurs.
//     */
//    public void write(String filePath) throws IOException {
//        PrintWriter writer = new PrintWriter(Files.newBufferedWriter(Paths.get(filePath)));
//        
//        // print dispatch solution
//        writer.printf("Dispatch solution\n");
//        for (int i = 0; i < problem.nTimeHorizons * problem.nRegions * (problem.nRegions - 1) / 2; i++) { 
//			writer.printf("%d ", dispatchSolution[i]);
//		}
//        writer.printf("\n\n");
//        
//        // print supply
//		int[] supply = new int[problem.nRegions];
//        for (int i = 0; i < problem.nRegions; i++) {
//        	supply[i] = problem.initialSupply[i];
//        }
//        for (int k = 0; k < problem.nTimeHorizons; k++) { 
//        	writer.printf("k = %d\n", k+1);
//        	// print dispatch solution
//        	/*
//        	writer.printf("Dispatch solution\n");
//			for (int i = 0; i < problem.nRegions; i++) {
//				for (int j = 0; j < problem.nRegions; j++) {
//					if (j < i) {
//						if (dispatchSolution[solutionMapping[k][j][i]] < 0) {
//							writer.printf("%d ", -dispatchSolution[solutionMapping[k][j][i]]);
//						} else {
//							writer.printf("0 ");
//						}
//					} else if (j == i) {
//						writer.printf("0 ");
//					} else if (j > i) {
//						if (dispatchSolution[solutionMapping[k][i][j]] > 0) {
//							writer.printf("%d ", dispatchSolution[solutionMapping[k][i][j]]);
//						} else {
//							writer.printf("0 ");
//						}
//					}
//				}
//				writer.printf("\n");
//			}
//			*/
//        	// print demand
//        	writer.printf("Demand   ");
//        	for (int i = 0; i < problem.nRegions; i++) {
//        		writer.printf("%d ", problem.predictedDemand[k][i]);
//        	}
//			// print supply after dispatch
//			for (int i = 0; i < problem.nRegions; i++) {
//				for (int j = 0; j < problem.nRegions; j++) {
//					if (i < j && dispatchSolution[solutionMapping[k][i][j]] > 0) {
//						supply[i] -= dispatchSolution[solutionMapping[k][i][j]];
//						supply[j] += dispatchSolution[solutionMapping[k][i][j]];
//					} else if (j < i && dispatchSolution[solutionMapping[k][j][i]] < 0) {
//						supply[i] -= -dispatchSolution[solutionMapping[k][j][i]];
//						supply[j] += -dispatchSolution[solutionMapping[k][j][i]];
//					}
//				}
//			}
//			writer.printf("\nSupply   ");
//        	for (int i = 0; i < problem.nRegions; i++) {
//        		writer.printf("%d ", supply[i]);
//        	}
//			// print supply after mobility
//			writer.printf("\nMobility ");
//			double[] temp = new double[problem.nRegions];
//			for (int ii = 0; ii < problem.nRegions; ii++) {
//				for (int jj = 0; jj < problem.nRegions; jj++) {
//					temp[ii] += supply[jj] * problem.mobilityPattern[jj][ii];
//				}
//			}
//			for (int ii = 0; ii < problem.nRegions; ii++) {
//				supply[ii] = (int) Math.round(temp[ii]);
//				writer.printf("%d ", supply[ii]);
//			}
//			writer.printf("\n\n");
//		}
//        
//        updateCost();
//        writer.printf("JD=%f\n", JDCost);
//        writer.printf("JE=%f\n", JECost);
//        writer.printf("Cost=%f\n", dispatchCost);
//        writer.close();
//    }
//}