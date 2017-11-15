package dsbp.model;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
//import java.nio.file.Files;
//import java.nio.file.Paths;

import dsbp.util.SimpleTokenizer;

/**
 * This class represents an Taxi demand and supply balancing problem.
 *
 * @author Hoang Tran
 */
public class Problem {

    /***
     * Number of regions
     */
    public final int nRegions;

    /***
     * Number of time horizons
     */
    public final int nTimeHorizons;

    /***
     * the weight factor of the objective function
     */
    public final float weightFactor;

    /***
     * Matrix with the distance of a region to region
     * distance[region][region]
     */
    public final double distances[][];

    /***
     * Matrix with the mobility pattern 
     * mobilityPattern[region][region]
     */
    public final double mobilityPattern[][];

    /***
     * Array with the initial supply 
     * initialSupply[region]
     */
    public final int initialSupply[];

    /***
     * Matrix with the predicted demand 
     * predictedDemand[time][region]
     */
    public final int predictedDemand[][];

    private int totalSupply;
    private int[] totalDemand;
    private double[] globalDemandSupplyRatio;
    
    /**
     * Instantiates a new Problem from a file.
     *
     * @param instancePath the instance file path
     */
    public Problem(String instancePath, int beta) throws IOException {
        BufferedReader reader =  new BufferedReader(new FileReader(instancePath)); //Files.newBufferedReader(Paths.get(instancePath));

        SimpleTokenizer token = new SimpleTokenizer(reader.readLine());

        // reading number of regions and number of time horizons
        nRegions = token.nextInteger();
        nTimeHorizons = token.nextInteger();
        weightFactor = beta;
        
        // initializing arrays
        distances = new double[nRegions][nRegions];
        mobilityPattern = new double[nRegions][nRegions];
        initialSupply = new int[nRegions];
        predictedDemand = new int[nTimeHorizons][nRegions];

        // skip next line
        reader.readLine();

        // reading distance
        for (int i = 0; i < nRegions; i++) {
            token = new SimpleTokenizer(reader.readLine());
            for (int j = 0; j < nRegions; j++) {
                distances[i][j] = token.nextFloat();
            }
        }

        // skip next line
        reader.readLine();

        // reading mobility pattern
        for (int i = 0; i < nRegions; i++) {
            token = new SimpleTokenizer(reader.readLine());
            for (int j = 0; j < nRegions; j++) {
                mobilityPattern[i][j] = token.nextFloat();
            }
        }
 
        // skip next line
        reader.readLine();

        // reading initial supply
        token = new SimpleTokenizer(reader.readLine());
        for (int i = 0; i < nRegions; i++) {
            initialSupply[i] = token.nextInteger();
        }

        // skip next line
        reader.readLine();

        // reading predicted demand
        for (int i = 0; i < nTimeHorizons; i++) {
            token = new SimpleTokenizer(reader.readLine());
            for (int j = 0; j < nRegions; j++) {
                predictedDemand[i][j] = token.nextInteger();
            }
        }

        // compute total supply
        totalSupply = 0;
        for (int i = 0; i < nRegions; i++) {
        	totalSupply += initialSupply[i];
        }
        
        // compute total demand
        totalDemand = new int[nTimeHorizons];
        for (int k = 0; k < nTimeHorizons; k++) {
        	totalDemand[k] = 0;
        	for (int i = 0; i < nRegions; i++) {
        		totalDemand[k] += predictedDemand[k][i];
        	}
        }
        
        // compute global demand supply ratio
        globalDemandSupplyRatio = new double[nTimeHorizons];
        for (int k = 0; k < nTimeHorizons; k++) {
        	globalDemandSupplyRatio[k] = totalDemand[k] * 1.0 / totalSupply;
        }

        reader.close();
    }

	public int getTotalDemand(int k) {
		
		return totalDemand[k];
	}

	public int getTotalSupply() {
		
		return totalSupply;
	}
	
	public double getGlobalDemandSupplyRatio(int k) {
		return globalDemandSupplyRatio[k];
	}
}
