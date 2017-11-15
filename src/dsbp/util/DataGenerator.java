package dsbp.util;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Random;

public class DataGenerator {
	
	private int nTimeHorizons;
	private int nRegions;
	@SuppressWarnings("unused")
	private int nRow;
	private int nCol;
	private int nSupply;
	private int nDemand;
	private double distances[][];
    private double mobilityPattern[][];
    private int initialSupply[];
    private int predictedDemand[][];
    
	public DataGenerator(int nTimeHorizons, int nRow, int nCol, int nSupply, int nDemand) {
		
		this.nTimeHorizons = nTimeHorizons;
		this.nRow = nRow;
		this.nCol = nCol;
		this.nRegions = nRow * nCol;
		this.nSupply = nSupply;
		this.nDemand = nDemand;
	}
	
	public void generate() {
		
		int[] randVec = new int[nRegions];
		int sum = 0;
		Random rand = new Random();
		
		// initializing arrays
        distances = new double[nRegions][nRegions];
        mobilityPattern = new double[nRegions][nRegions];
        initialSupply = new int[nRegions];
        predictedDemand = new int[nTimeHorizons][nRegions];
        
        // generate distances
        int x1 = 0;
        int y1 = 0;
        int x2 = 0;
        int y2 = 0;
        for (int i = 0; i < nRegions; i++) {
        	for (int j = 0; j < nRegions; j++) {
        		x1 = i / nCol;
        		y1 = i % nCol;
        		x2 = j / nCol;
        		y2 = j % nCol;
        		distances[i][j] = Math.round(Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) * 100) / 100.0;
        	}
        }
        
        // generate mobility pattern
        for (int i = 0; i < nRegions; i++) {
        	sum = 0;
        	for (int j = 0; j < nRegions; j++) {
        		randVec[j] = rand.nextInt(100);
        		sum += randVec[j];
        	}
        	for (int j = 0; j < nRegions; j++) {
        		mobilityPattern[i][j] = randVec[j] * 1.0 / sum;
        	}
        }
        
        // generate initial supply
        sum = 0;
        for (int i = 0; i < nRegions; i++) {
        	randVec[i] = rand.nextInt(100);
    		sum += randVec[i];
        }
        for (int i = 0; i < nRegions; i++) {
        	initialSupply[i] = (int) Math.round(randVec[i] * 1.0 / sum * nSupply);
        }
        
        // generate predicted demand
        for (int i = 0; i < nTimeHorizons; i++) {
        	sum = 0;
            for (int j = 0; j < nRegions; j++) {
            	randVec[j] = rand.nextInt(100);
        		sum += randVec[j];
            }
            for (int j = 0; j < nRegions; j++) {
                predictedDemand[i][j] = (int) Math.round(randVec[j] * 1.0 / sum * nDemand);
            }
        }
	}
	
	/**
     * Writes the solution to a file.
     *
     * @param filePath the output file path.
     * @throws IOException in case any IO error occurs.
     */
    public void write(String filePath) throws IOException {
        PrintWriter writer = new PrintWriter(Files.newBufferedWriter(Paths.get(filePath)));

        // writing number of regions and number of time horizons
        writer.printf("%d %d", nRegions, nTimeHorizons);
        writer.printf("\n\n");
        
        // writing distances
        for (int i = 0; i < nRegions; i++) {
        	for (int j = 0; j < nRegions; j++) {
        		writer.printf("%f ", distances[i][j]);
        	}
        	writer.printf("\n");
        }
        writer.printf("\n");
        
        // writing mobility pattern
        for (int i = 0; i < nRegions; i++) {
        	for (int j = 0; j < nRegions; j++) {
        		writer.printf("%f ", mobilityPattern[i][j]);
        	}
        	writer.printf("\n");
        }
        writer.printf("\n");
        
        // writing initial supply
        for (int i = 0; i < nRegions; i++) {
        	writer.printf("%d ", initialSupply[i]);
        }
        writer.printf("\n\n");
        
        // writing predicted demand
        for (int i = 0; i < nTimeHorizons; i++) {
        	for (int j = 0; j < nRegions; j++) {
        		writer.printf("%d ", predictedDemand[i][j]);
        	}
        	writer.printf("\n");
        }

        writer.close();
    }
}
