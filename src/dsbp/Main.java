package dsbp;

import java.io.*;
import java.util.*;

import dsbp.algorithm.initialization.*;
import dsbp.algorithm.heuristic.*;
import dsbp.model.*;
import dsbp.util.*;

/**
 * This class is the Main class of the program, responsible of parsing the input, instantiating moves and heuristics and
 * printing the results.
 *
 * @author Hoang Tran
 */
public class Main {

    public static long startTimeMillis = System.currentTimeMillis();

    public static String algorithm = "esa";
    public static String inFile;
    public static String outFile;
    public static long seed = 0;
    public static long maxIters = ( long ) 1e8;
    public static long timeLimit = 3600 * 1000;
    
    public static void main(String[] args) throws IOException {
        Locale.setDefault(new Locale("en-US"));
        readArgs(args);
        seed = System.currentTimeMillis();
        Problem problem = new Problem(inFile, 100);
        Random random = new Random(seed);

        Heuristic solver = null;
        switch (algorithm) {
            case "esa":
                solver = new ESA(problem, random);
                break;
            case "descent":
            	solver = new Descent(problem, random);
            	break;
            default:
                System.exit(-1);
                return;
        }

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
        assert solution.validate();
        solution = solver.run(solution, timeLimit, maxIters, System.out);
        solution.validate();

        System.out.printf("    \\--------------------------------------------------------/\n\n");

        System.out.printf("Best makespan.....: %f\n", solution.getTotalCost());
        System.out.printf("N. of Iterations..: %d\n", solver.getNIters());
        System.out.printf("Total runtime.....: %.2fs\n", (System.currentTimeMillis() - startTimeMillis) / 1000.0);
        solution.updateTotalCost();
        System.out.println(solution.getTotalCost());
        solution.write(outFile);
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
        System.out.println("    -seed <seed>           : random seed (default: " + seed + ").");
        System.out.println("    -maxIters <maxIters>   : maximum number of consecutive rejections (default: Long.MAXVALUE).");
        System.out.println("    -time <timeLimit>      : time limit in seconds (default: " + timeLimit + ").");
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
                default:
                    printUsage();
                    System.exit(-1);
            }
        }
    }
}
