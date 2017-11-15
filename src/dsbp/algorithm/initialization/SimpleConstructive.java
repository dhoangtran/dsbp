package dsbp.algorithm.initialization;

import java.util.*;

import dsbp.model.*;

/**
 * This class contains simple constructive procedures for the UPMSP.
 *
 * @author DuyHoang Tran
 */
public class SimpleConstructive {

    /**
     *
     * @param problem problem
     * @param random  random number generator
     * @return the solution generated
     */
    public static Solution randomSolution(Problem problem, Random random) {
        Solution solution = new Solution(problem);
        
        for (int i = 0; i < solution.size; i++) {
        	solution.dispatch[i] = 0;
        }

        solution.updateTotalCost();
        return solution;
    }
}
