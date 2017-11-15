package dsbp.util;

import dsbp.*;
import dsbp.model.*;
import dsbp.algorithm.neighborhood.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

/**
 * Simple class with some util methods.
 *
 * @author Tulio Toffolo
 */
public class Util {

    /**
     * Convert a Long number to a readable string.
     *
     * @param value long number.
     * @return the "human-friendly" (readable) number as a String.
     */
    public static String longToString(long value) {
        return value >= 1e11 ? String.format("%.0fG", value / 1e9)
          : value >= 1e7 ? String.format("%.0fM", value / 1e6)
          : (value >= 1e4) ? String.format("%.0fK", value / 1e3)
          : Long.toString(value);
    }

    /**
     * Calls printf after checking if the PrintStream is not null.
     *
     * @param output the output stream.
     * @param format the String with the format.
     * @param args   "printf" arguments.
     */
    public static void safePrintf(PrintStream output, String format, Object... args) {
        if (output != null) {
            output.printf(format, args);
        }
    }

    /**
     * Prints the statistics of a Move using the table stable, after checking that the PrintStream is not null.
     *
     * @param output  the output stream.
     * @param move    the Move considered.
     * @param special some informative String to print after the row.
     */
    public static void safePrintMoveStatistics(PrintStream output, Move move, String special) {
        if (output != null) {
            output.printf("    | %-18s | %8s | %8s | %8s | %8s |\n",
              move.name,
              longToString(move.getNImprovements()),
              longToString(move.getNSideways()),
              longToString(move.getNAccepts()),
              longToString(move.getNRejects())
            );
        }
    }


    /**
     * Prints the current solution status after checking that the PrintStream is not null.
     *
     * @param output       the output stream.
     * @param nIters       the current iteration number.
     * @param bestSolution the best solution object.
     * @param solution     the current solution object.
     * @param special      some informative String to print after the row.
     */
    public static void safePrintStatus(PrintStream output, long nIters, Solution bestSolution, Solution solution, String special) {
        if (output != null) {
            output.printf("    | %8s | %8s | %8f | %8f | %10.2f | %s\n",
              longToString(nIters),
              Main.bestKnown == Integer.MAX_VALUE ? "-" : String.format("%8.2f",
                100. * ( double ) (solution.getTotalCost() - Main.bestKnown) / ( double ) Main.bestKnown),
              bestSolution.getTotalCost(), solution.getTotalCost(),
              (System.currentTimeMillis() - Main.startTimeMillis) / 1000.0, special
            );
        }
    }

    /**
     * Prints the text maintaining the table style,after checking that the PrintStream is not null.
     *
     * @param output  the output stream.
     * @param text    text to print inside the table.
     * @param special some informative String to print after the row.
     */
    public static void safePrintText(PrintStream output, String text, String special) {
        if (output != null) {
            output.printf("    | %-40s | %10.2f | %s\n", text, (System.currentTimeMillis() - Main.startTimeMillis) / 1000.0, special);
        }
    }
    
    public static void writeToFile(String outFile, String text) throws IOException {
    	PrintWriter writer;
		writer = new PrintWriter(Files.newBufferedWriter(Paths.get(outFile), StandardCharsets.US_ASCII, StandardOpenOption.CREATE, StandardOpenOption.APPEND));
		writer.printf(text);
		writer.close();
    }
}
