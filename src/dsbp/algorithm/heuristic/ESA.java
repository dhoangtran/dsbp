package dsbp.algorithm.heuristic;

import dsbp.algorithm.neighborhood.ESAMove;
import dsbp.model.*;
import dsbp.util.*;

import java.io.*;
import java.util.*;

/**
 * This class is a Simulated Annealing implementation.
 *
 * @author Hoang Tran
 */
public class ESA extends Heuristic {

    /**
     * SA parameters.
     */
    //private double alpha, t0;
    //private int saMax = 10000;
    //private final static double EPS = 0.01;//1e-6;
	private final static double RATMAX = 0.02; 
	private final static double RATMIN = 0.005; 
	private final static double EXTSTP = 1.2; 
	private final static double SHRSTP = 0.8;
	// lower and upper bounds
	private int[] X0MIN;
	private int[] X0MAX;
	// end of temperature state parameters 
	private int N1;
	private int N2;
	private double TMPINI;
	// stopping tests parameters
	private double EPSREL = 1.0D-06;
	private double EPSABS = 1.0D-06;//1.0D-08;
	private double NFMAX = 5000;
	// temperature adjustment rule parameters
	private final static double RMXTMP = 0.99;
	private final static double RMITMP = 0.1;
	// step vector adjustment rule parameters
	private int[] STPMST;
	private double[] ROSTEP;
	//private int RATMIN;
	//private int RATMAX;
	//private int EXTSTP;
	//private int SHRSTP;
	// number of FOBJ evaluations
	private int NFOBJ;
	private int MVOKST;
	private int[] MOKST;
	private int NMVUST;
	private int NMVST;
	private double SDGYUP;
	private double ELOWST;
	private int[] MTOTST;
	private double AVGYST;
    /**
     * Instantiates a new SA.
     *
     * @param problem problem reference
     * @param random  random number generator.
     * @param alpha   cooling rate for the simulated annealing
     * @param t0      initial temperature, T0
     * @param saMax   number of iterations before update the temperature
     */
    public ESA(Problem problem, Random random) {
        super(problem, random, "SA");

        // initializing simulated annealing parameters
        
    }

    /**
     * Executes the Simulated Annealing.
     *
     * @param initialSolution the initial (input) solution.
     * @param timeLimitMillis the time limit (in milliseconds).
     * @param maxIters        the maximum number of iterations without improvements to execute.
     * @param output          output PrintStream for logging purposes.
     * @return the best solution encountered by the SA.
     */
    public Solution run(Solution initialSolution, long timeLimitMillis, long maxIters, PrintStream output) {
    	ESAMove move = new ESAMove(problem, random, 1);
    	long finalTimeMillis = System.currentTimeMillis() + timeLimitMillis;

        bestSolution = initialSolution;
        Solution solution = initialSolution.clone();

        // Step 1: Initializations
        N2 = 150 * solution.size;
        
        double avgRange = 0;
        for (int i = 0; i < solution.size; i++) {
        	avgRange += solution.getRange(i);	
        }
        avgRange /= solution.size;
        int step = (int) Math.round(avgRange * 0.10);
        STPMST = new int[solution.size];
        Arrays.fill(STPMST, step);
        
        int count = 0;
        ArrayList<Double> deltas = new ArrayList<Double>();
        Solution init = solution.clone();
        for (int i = 0; i < 1000; i++) {
        	double delta = move.doMove(init, STPMST);
        	if (delta < 0) {
        		acceptMove(move);
        	} else {
        		count++;
        		deltas.add(delta);
        		rejectMove(move);
        	}
        }
        Collections.sort(deltas);
        double avg = 0;
        for (int i = 0; i < 10; i++) {
        	avg += deltas.get(i);
        }
        avg = avg / 10;
        TMPINI = -avg/Math.log(0.5);
        double temperature = TMPINI;
        
        NFOBJ = 1;
        MVOKST = 0;
        NMVUST = 0;
        NMVST = 0;
        MOKST = new int[solution.size];
        MTOTST = new int[solution.size];
        ELOWST = solution.getTotalCost();
        AVGYST = 0;
        SDGYUP = 0;
        
        int nItersWithoutImprovement = 0;
        maxIters = 1500 * solution.size;
        while (System.currentTimeMillis() < finalTimeMillis && nItersWithoutImprovement++ < maxIters) {
            double delta = move.doMove(solution, STPMST);
            MTOTST[move.position]++;
            AVGYST = AVGYST + solution.getTotalCost();
            NMVST++;
            NFOBJ++;
            
            // if solution is improved...
            if (delta < 0) {
                acceptMove(move);
                MOKST[move.position]++;
            	if (solution.getTotalCost() < bestSolution.getTotalCost()) {
            		nItersWithoutImprovement = 0;
                    bestSolution = solution.clone();
                    Util.safePrintStatus(output, nIters, bestSolution, solution, "*");
                }
            	if (solution.getTotalCost() < ELOWST) ELOWST = solution.getTotalCost();
            }

            // solution is not improved, but may be accepted with a probability...
            else {
                double x = random.nextDouble();
                if (x < 1 / Math.exp(delta / temperature)) {
                    acceptMove(move);
                    MOKST[move.position]++;
                    NMVUST++;
                }

                // if solution is rejected..
                else {
                    rejectMove(move);
                }
                
                SDGYUP = SDGYUP + delta;
            }

            // if necessary, updates temperature
            if (NMVST > N2) {
            	AVGYST = AVGYST/NMVST;
            	double RFTMP = Math.max(Math.min(ELOWST/AVGYST, RMXTMP), RMITMP);
                temperature = RFTMP * temperature;
            	System.out.println(temperature);
                // Step 7: Step Vector Adjustment
                for (int i = 0; i < solution.size; i++) {
                	double ROK = MOKST[i] * 1.0 / MTOTST[i];
                	if (ROK > RATMAX) STPMST[i] = (int) Math.round(STPMST[i] * EXTSTP);
                	else if (ROK < RATMIN) STPMST[i] = (int) Math.round(STPMST[i] * SHRSTP);
                }
                
                // Step 9: Initialization of a New Temperature Stage
                MVOKST = 0;
                NMVUST = 0;
            	NMVST = 0;
            	Arrays.fill(MOKST, 0);
            	Arrays.fill(MTOTST, 0);
            	ELOWST = solution.getTotalCost();
            	AVGYST = 0;
            	SDGYUP = 0;
            }
            
            nIters++;
        }

        return bestSolution;
    }

    /**
     * Returns the string representation of this heuristic.
     *
     * @return the string representation of this heuristic (with parameters values).
     */
    public String toString() {
        return String.format("Simulated Annealing");
    }


    //private void estimateT0(Solution initialSolution, int nNeighbors, double ratio) {
        //Solution solution = initialSolution.clone();
        //List<Integer> neighborValues = new ArrayList<>(nNeighbors);
        //
        //for (int i = 0; i < nNeighbors; i++) {
        //    Move move = selectMove(solution);
        //    int delta = move.doMove(solution);
        //
        //    neighborValues.add(delta);
        //}
        //neighborValues.sort(Integer::compare);
        //
        //int t = 1;
        //int
        //
        //int idealDelta = neighborValues.get(( int ) (nNeighbors * ratio));
        //1/FastMath.log(delta / temperature);
        //
        //1/ratio
    //}
}
