package dsbp.algorithm.heuristic;

import java.io.*;
import java.util.*;

import dsbp.algorithm.neighborhood.*;
import dsbp.model.*;

/**
 * This abstract class represents a Heuristic (or Local Search method). The basic methods and neighborhood selection are
 * included.
 *
 * @author Hoang Tran
 */
public abstract class Heuristic {

    public final Problem problem;
    public final Random random;
    public final String name;

    protected Solution bestSolution;
    protected long nIters = 0;

    /**
     * Instantiates a new Heuristic.
     *
     * @param problem the problem reference.
     * @param random  the random number generator.
     * @param name    the name
     */
    public Heuristic(Problem problem, Random random, String name) {
        this.problem = problem;
        this.random = random;
        this.name = name;
    }

    /**
     * Accepts move and updates learning algorithm (if present).
     *
     * @param move the move to be accepted.
     */
    public void acceptMove(Move move) {
        move.accept();
    }

    /**
     * Rejects move and updates learning algorithm (if present).
     *
     * @param move the move to be rejected.
     */
    public void rejectMove(Move move) {
        move.reject();
    }

    /**
     * Runs the local search, returning the best solution obtained..
     *
     * @param solution        the initial (input) solution.
     * @param timeLimitMillis the time limit in milliseconds.
     * @param maxIters        the maximum number of iterations to execute.
     * @param output          the output
     * @return the solution
     */
    public abstract Solution run(Solution solution, long timeLimitMillis, long maxIters, PrintStream output);

    /**
     * Gets best solution.
     *
     * @return the best solution obtained so far.
     */
    public Solution getBestSolution() {
        return bestSolution;
    }

    /**
     * Gets the number of iterations executed.
     *
     * @return the n iters
     */
    public long getNIters() {
        return nIters;
    }

    /**
     * Returns the string representation of the heuristic.
     *
     * @return the string representation of the heuristic.
     */
    public String toString() {
        return name;
    }
}
