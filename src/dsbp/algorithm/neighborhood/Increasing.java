package dsbp.algorithm.neighborhood;

import java.util.*;

import dsbp.model.*;

/**
 * This class represents a Shift Move (as in the paper). A neighbor in the Shift
 * Move is generated by re-scheduling one job from the machine with the largest
 * total execution time (or a random machine) to another position in the
 * machine. The parameter "useMakespanMachine" determines whether the machine
 * with the largest total execution time is always used.
 *
 * @author Hoang Tran
 */
public class Increasing extends Move {

	private int pos;
	private int nSupply;
    /**
     * Instantiates a new Shift Move.
     *
     * @param problem            problem.
     * @param random             random number generator.
     * @param priority           the priority of this neighborhood.
     * @param useMakespanMachine true if the makespan machine should be always
     *                           considered or false otherwise.
     */
    public Increasing(Problem problem, Random random, int priority) {
        super(problem, random, "Increasing", priority);
    }

    public void accept() {
        super.accept();
    }

    public double doMove(Solution solution) {
        super.doMove(solution);

        boolean valid = false;
        do {
        	pos = random.nextInt(solution.size);
        	//nSupply = random.nextInt(3) + 1;
        	nSupply = 1;
        	solution.dispatch[pos] += nSupply;
        	if (solution.validate() == true) {
        		valid = true;
        	} else {
        		solution.dispatch[pos] -= nSupply;
        	}
        } while (!valid);
       
        solution.updateTotalCost();
        return deltaCost = solution.getTotalCost() - initialCost;
    }

    public boolean hasMove(Solution solution) {
        return solution.validate();
    }

    public void reject() {
        super.reject();
        
        currentSolution.dispatch[pos] -= nSupply;
        currentSolution.updateTotalCost();
    }
}