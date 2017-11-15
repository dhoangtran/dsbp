package dsbp;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;

import dsbp.algorithm.initialization.*;
import dsbp.algorithm.heuristic.*;
import dsbp.algorithm.neighborhood.*;
import dsbp.model.*;
import dsbp.util.*;

/**
 * This class is the Main class of the program, responsible of parsing the input, instantiating moves and heuristics and
 * printing the results.
 *
 * @author DuyHoang Tran
 */
public class Main {

    // region solver parameters and default values

    public static long startTimeMillis = System.currentTimeMillis();

    public static String algorithm = "esa";
    public static String inFile;
    public static String outFile = null;

    public static long seed = 0;
    public static long maxIters = ( long ) 10000000;//1e8;
    public static long timeLimit = 3600 * 1000;

    public static int bestKnown = Integer.MAX_VALUE;

    // ILS
    public static long rnaMax = 100000;
    public static int itersP = 20;
    public static int p0 = 4;//80;
    public static int pMax = 6;

    // LAHC
    public static int listSize = ( int ) 1000;

    // SA (Simulated Annealing)
    public static double alpha = 0.95;//0.99;
    public static int saMax = ( int ) 100000;//1e7;
    public static double t0 = 0.5;//1;

    // SCHC
    public static int stepSize = 1000;
    
    // TS
    public static int tabuListSize = 10;
    
    // Neighborhoods
    public static boolean neighborhoods[];
    
    static {
        neighborhoods = new boolean[6 * 4];
        for (int i = 0; i < neighborhoods.length; i++)
            neighborhoods[i] = true;
    }

    // endregion solver parameters and default values

    /**
     * The entry point of application.
     *
     * @param args the input arguments
     * @throws IOException if any IO error occurs.
     */
    public static void main(String[] args) throws IOException {
        Locale.setDefault(new Locale("en-US"));
        if (args.length < 2) {
            System.out.printf("Need algorithm parameter.");
            System.exit(-1);
            return;
        }
        
        String dataPath = args[0];
        algorithm = args[1];
        DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd");
        Date date = new Date();
        outFile = String.format("%s/out/%s_%s.txt", dataPath, algorithm, dateFormat.format(date));
        String[] sizes = {"16x4", "25x4", "36x4", "50x4"};
        //String[] sizes = {"50x4"};
        for (String size : sizes) {    		
    		int[] betas = {10, 20, 50, 100};
		    for (int beta : betas) {
		    	//Util.writeToFile(outFile, String.format("size=%s,beta=%d\n", size, beta));
		    	Util.writeToFile(outFile, "\n");
		    	for (int id = 1; id <= 10; id++) {
		    		System.out.printf("algorithm=%s,size=%s,beta=%d,id=%d\n", algorithm, size, beta, id);
		    		inFile = String.format("%s/in/%s_%d.txt", dataPath, size, id);
		    		Problem problem = new Problem(inFile, beta);
		        	Solution bestSolution = new Solution(problem);
		        	for (int k = 0; k < 1; k++) {
		        		seed = System.currentTimeMillis();
			        	Random random = new Random(seed);
		        		Heuristic solver = null;
		        		switch (algorithm) {
		                	case "ils":
		                		solver = new ILS(problem, random, rnaMax, itersP, p0, pMax);
		                		break;
		                	case "lahc":
		                		solver = new LAHC(problem, random, listSize);
		                		break;
		                	case "schc":
		                		solver = new SCHC(problem, random, stepSize);
		                		break;
		                	case "sa":
		                		solver = new SA(problem, random, alpha, t0, saMax);
		                		break;
		                	case "esa":
		                        solver = new ESA(problem, random);
		                        break;
		                	case "ts":
		                		//solver = new Tabu(problem, random, tabuListSize);
		                		//break;
		                	default:
		                        System.exit(-1);
		                        return;
		        		}    
			        	// adding moves (neighborhoods)
			        	//createNeighborhoods(problem, random, solver);
			
			        	// re-starting time counting (after reading files)
			        	startTimeMillis = System.currentTimeMillis();
			
			        	// generating initial solution
			        	Solution solution = new Solution(problem); //SimpleConstructive.randomSolution(problem, random);
			        	//assert solution.validate();
			
			        	// running stochastic local search
			        	//if (solver.getMoves().size() > 0)
			        	solution = solver.run(solution, timeLimit, maxIters, System.out);
			        	//assert solution.validate();
			        	solution.updateTotalCost();
			        	if (solution.getTotalCost() < bestSolution.getTotalCost()) bestSolution = solution;
		        	}
		        	Util.writeToFile(outFile, String.format("%d %d %d %d\n", 
		        			Math.round(bestSolution.getDrivingCost()), 
		        			Math.round(bestSolution.getBalancingCost()), 
		        			Math.round(bestSolution.getTotalCost()), 
		        			(System.currentTimeMillis() - startTimeMillis) / 1000));
		        }
        	}
        }
    }

    public static void main2(String[] args) throws IOException {
        Locale.setDefault(new Locale("en-US"));
        algorithm = "ils";
        String size = "16x4";
        int id = 1;
        int beta = 10;
    	
        inFile = String.format("/Users/tdhoang/Dropbox/Research/Taxi dispatching/Datasets/in/%s_%d.txt", size, id);
    	outFile = "/Users/tdhoang/Dropbox/Research/Taxi dispatching/Datasets/out/test.txt";
    	Util.writeToFile(outFile, String.format("%s, %s, $d\n", algorithm, size, beta));
    	Problem problem = new Problem(inFile, beta);
    	Solution bestSolution = new Solution(problem);
    	Util.writeToFile(outFile, String.format("%d %d %d\n", Math.round(bestSolution.getDrivingCost()), 
    			Math.round(bestSolution.getBalancingCost()), Math.round(bestSolution.getTotalCost())));
    	for (int k = 0; k < 1; k++) {
    		seed = System.currentTimeMillis();
    		Random random = new Random(seed);
    		Heuristic solver = null;
    		// ILS parameter
    	    
    		
    		//solver = new ILS(problem, random, rnaMax, itersP, p0, pMax);
    		//solver = new LAHC(problem, random, listSize);
    		//solver = new SCHC(problem, random, stepSize);
    		//solver = new SA(problem, random, alpha, t0, saMax);
        	//solver = new Descent(problem, random);        
        	//solver = new Tabu(problem, random, 10);
                
        	// adding moves (neighborhoods)
        	//createNeighborhoods(problem, random, solver);

        	// re-starting time counting (after reading files)
        	startTimeMillis = System.currentTimeMillis();

        	// generating initial solution
        	Solution solution = new Solution(problem);//SimpleConstructive.randomSolution(problem, random);
        	//assert solution.validate();

        	// running stochastic local search
        	//if (solver.getMoves().size() > 0)
        	//solution = solver.run(solution, timeLimit, maxIters, System.out);
        	solution = solver.run(solution, timeLimit, maxIters, System.out);
        	assert solution.validate();
        	solution.updateTotalCost();
        	if (solution.getTotalCost() < bestSolution.getTotalCost()) bestSolution = solution;
        	Util.writeToFile(outFile, String.format("%d %d %d\n", Math.round(solution.getDrivingCost()), Math.round(solution.getBalancingCost()), Math.round(solution.getTotalCost())));
    	}
    	Util.writeToFile(outFile, String.format("best: %d %d %d\n", Math.round(bestSolution.getDrivingCost()), Math.round(bestSolution.getBalancingCost()), Math.round(bestSolution.getTotalCost())));
    }

    public static void main1(String[] args) throws IOException {
        Locale.setDefault(new Locale("en-US"));
        readArgs(args);
        seed = System.currentTimeMillis();
        Problem problem = new Problem(inFile, 100);
        Random random = new Random(seed);

        Heuristic solver = null;
        switch (algorithm) {
            case "lahc":
                solver = new LAHC(problem, random, listSize);
                break;
            case "lahc-ils":
                solver = new ILS(problem, random, new LAHC(problem, random, listSize), rnaMax, itersP, p0, pMax);
                break;
            case "ils":
                solver = new ILS(problem, random, rnaMax, itersP, p0, pMax);
                break;
            case "sa":
                solver = new SA(problem, random, alpha, t0, saMax);
                break;
            case "esa":
                solver = new ESA(problem, random);
                break;
            case "sa-ils":
                solver = new ILS(problem, random, new SA(problem, random, alpha, t0, saMax), rnaMax, itersP, p0, pMax);
                break;
            case "schc":
                solver = new SCHC(problem, random, stepSize);
                break;
            case "schc-ils":
                solver = new ILS(problem, random, new SCHC(problem, random, stepSize), rnaMax, itersP, p0, pMax);
                break;
            case "descent":
            	solver = new Descent(problem, random);
            	break;
            case "ts":
            	//solver = new Tabu(problem, random, 10);
            	//break;
            case "es":
            	//solver = new ES(problem, random);
            	//break;
            	
            default:
                System.exit(-1);
                return;
        }

        // adding moves (neighborhoods)
        createNeighborhoods(problem, random, solver);

        System.out.printf("Instance....: %s\n", inFile);
        System.out.printf("Algorithm...: %s\n", solver);
        System.out.printf("Other params: maxIters=%s, seed=%d, timeLimit=%.2fs\n\n", Util.longToString(maxIters), seed, timeLimit / 1000.0);
        System.out.printf("    /--------------------------------------------------------\\\n");
        System.out.printf("    | %8s | %8s | %8s | %8s | %10s | %s\n", "Iter", "RDP(%)", "S*", "S'", "Time", "");
        System.out.printf("    |----------|----------|----------|----------|------------|\n");

        // re-starting time counting (after reading files)
        startTimeMillis = System.currentTimeMillis();

        // generating initial solution
        Solution solution = SimpleConstructive.randomSolution(problem, random);
        //Solution solution = new Solution(problem);
		//solution.readFromFile("C:/Users/r0660215/workspace/balancing/data/out/50x4_1_out_cvx.txt");
        Util.safePrintStatus(System.out, 0, solution, solution, "s0");
        assert solution.validate();

        // running stochastic local search
        if (solver.getMoves().size() > 0)
        	//maxIters = 1000 * solution.size;
            solution = solver.run(solution, timeLimit, maxIters, System.out);
        solution.validate();

        System.out.printf("    \\--------------------------------------------------------/\n\n");

        System.out.printf("Neighborhoods statistics (values in %%):\n\n");
        System.out.printf("    /----------------------------------------------------------------\\\n");
        System.out.printf("    | %-18s | %8s | %8s | %8s | %8s |\n", "Move", "Improvs.", "Sideways", "Accepts", "Rejects");
        System.out.printf("    |--------------------|----------|----------|----------|----------|\n");
        for (Move move : solver.getMoves())
            Util.safePrintMoveStatistics(System.out, move, "");
        System.out.printf("    \\----------------------------------------------------------------/\n\n");

        if (bestKnown != Integer.MAX_VALUE)
            System.out.printf("Best RDP..........: %.4f%%\n", 100 * ( double ) (solution.getTotalCost() - bestKnown) / ( double ) bestKnown);
        System.out.printf("Best makespan.....: %f\n", solution.getTotalCost());
        System.out.printf("N. of Iterations..: %d\n", solver.getNIters());
        System.out.printf("Total runtime.....: %.2fs\n", (System.currentTimeMillis() - startTimeMillis) / 1000.0);
        solution.updateTotalCost();
        System.out.println(solution.getTotalCost());
        solution.write(outFile);
    }

    private static void createNeighborhoods(Problem problem, Random random, Heuristic solver) {
        int index = -1;

        if (neighborhoods[++index]) solver.addMove(new Increasing(problem, random, 1));
        if (neighborhoods[++index]) solver.addMove(new Descreasing(problem, random, 1));
    }

    /**
     * Prints the program usage.
     */
    public static void printUsage() {
        System.out.println("Usage: java -jar upmsp.jar <input> <output> [options]");
        System.out.println("    <input>  : Path of the problem input file.");
        System.out.println("    <output> : Path of the (output) solution file.");
        System.out.println();
        System.out.println("Options:");
        System.out.println("    -algorithm <algorithm> : ils, lahc, lahc-ils, sa, sa-ils, schc or schc-ils (default: " + algorithm + ").");
        System.out.println("    -bestKnown <makespan>  : best known makespan for RDP output (default: " + bestKnown + ").");
        System.out.println("    -seed <seed>           : random seed (default: " + seed + ").");
        System.out.println("    -maxIters <maxIters>   : maximum number of consecutive rejections (default: Long.MAXVALUE).");
        System.out.println("    -time <timeLimit>      : time limit in seconds (default: " + timeLimit + ").");
        System.out.println();
        System.out.println("    ILS parameters:");
        System.out.println("        -rnamax <rnamax> : maximum rejected iterations in the descent phase of ILS (default: " + rnaMax + ").");
        System.out.println("        -itersP <itersP> : number of iterations per perturbation level for ILS (default: " + itersP + ").");
        System.out.println("        -p0 <p0>         : initial perturbation level for ILS (default: " + p0 + ").");
        System.out.println("        -pMax <pMax>     : maximum steps up (each step of value p0) for ILS perturbation's level (default: " + pMax + ").");
        System.out.println();
        System.out.println("    LAHC parameters:");
        System.out.println("        -listSize <listSize> : LAHC list size  (default: " + listSize + ").");
        System.out.println();
        System.out.println("    SA parameters:");
        System.out.println("        -alpha <alpha> : cooling rate for the Simulated Annealing (default: " + alpha + ").");
        System.out.println("        -samax <samax> : iterations before updating the temperature for Simulated Annealing (default: " + saMax + ").");
        System.out.println("        -t0 <t0>       : initial temperature for the Simulated Annealing (default: " + t0 + ").");
        System.out.println();
        System.out.println("    SCHC parameters:");
        System.out.println("        -stepSize <stepSize> : SCHC step size (default: " + stepSize + ").");
        System.out.println();
        System.out.println("    Neighborhoods selection:");
        System.out.println("        -n <id,policy,value> : disables a policy(0..3) for neighborhood id(0..5) if value = 0 and enables it otherwise.");
        System.out.println();
        System.out.println("Examples:");
        System.out.println("    java -jar upmsp.jar instance.txt solution.txt");
        System.out.println("    java -jar upmsp.jar instance.txt solution.txt -algorithm sa -alpha 0.98 -samax 1000 -t0 100000");
        System.out.println();
    }

    /**
     * Reads the input arguments.
     *
     * @param args the input arguments
     */
    public static void readArgs(String args[]) {
        if (args.length < 2) {
            printUsage();
            System.exit(-1);
        }

        int index = -1;

        inFile = args[++index];
        outFile = args[++index];

        while (index < args.length - 1) {
            String option = args[++index].toLowerCase();

            switch (option) {
                case "-algorithm":
                    algorithm = args[++index].toLowerCase();
                    break;
                case "-seed":
                    seed = Integer.parseInt(args[++index]);
                    break;
                case "-maxiters":
                    maxIters = Long.parseLong(args[++index]);
                    break;
                case "-time":
                    timeLimit = Math.round(Double.parseDouble(args[++index]) * 1000.0);
                    break;

                case "-bestknown":
                    bestKnown = Integer.parseInt(args[++index]);
                    break;

                // ILS
                case "-rnamax":
                    rnaMax = Long.parseLong(args[++index]);
                    break;
                case "-itersp":
                    itersP = Integer.parseInt(args[++index]);
                    break;
                case "-p0":
                    p0 = Integer.parseInt(args[++index]);
                    break;
                case "-pmax":
                    pMax = Integer.parseInt(args[++index]);
                    break;

                // LAHC
                case "-listsize":
                    listSize = Integer.parseInt(args[++index]);
                    break;

                // SA
                case "-alpha":
                    alpha = Double.parseDouble(args[++index]);
                    break;
                case "-samax":
                    saMax = Integer.parseInt(args[++index]);
                    break;
                case "-t0":
                    t0 = Double.parseDouble(args[++index]);
                    break;

                // SCHC
                case "-stepsize":
                    stepSize = Integer.parseInt(args[++index]);
                    break;

                // Neighborhoods selection
                case "-n":
                    String[] values = args[++index].split(",");
                    int i = Integer.parseInt(values[0]) * 4 + Integer.parseInt(values[1]);
                    neighborhoods[i] = values[2].equals("1");
                    break;

                default:
                    printUsage();
                    System.exit(-1);
            }
        }
    }
}
