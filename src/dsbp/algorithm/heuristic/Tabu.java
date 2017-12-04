package dsbp.algorithm.heuristic;


/**
 * Basic fixed size implementation for Tabu list
 * @author Hoang Tran
 *
 */
//class TabuList implements Iterable<Solution> {
//	
//	private CircularFifoQueue<Solution> tabuList;
//	
//	/**
//	 * Construct a new {@link TabuList}
//	 * @param size the size of the tabu list
//	 */
//	public TabuList(Integer size) {
//		this.tabuList = new CircularFifoQueue<Solution>(size);
//	}
//
//	public void add(Solution solution) {
//		tabuList.add(solution);
//	}
//	
//	public int size() {
//		return tabuList.size();
//	}
//
//	public Boolean contains(Solution solution) {
//		return tabuList.contains(solution);
//	}
//	
//	@Override
//	public Iterator<Solution> iterator() {
//		return tabuList.iterator();
//	}
//
//	/**
//	 * This method does not perform any update in the tabu list,
//	 * due to the fixed size nature of this implementation
//	 */
//	public void updateSize(Integer currentIteration, Solution bestSolutionFound) {
//		//Do nothing, this implementation has a fixed size
//	}
//}
//
///**
// * This class is a Tabu Search implementation.
// *
// * @author Hoang Tran
// */
//public class Tabu extends Heuristic {
//
//    /**
//     * Tabu parameters.
//     */
//	private int tabuListSize;
//	private TabuList tabuList;
//	/**
//     * Instantiates a new SA.
//     *
//     * @param problem problem reference
//     * @param random  random number generator.
//     * @param alpha   cooling rate for the simulated annealing
//     * @param t0      initial temperature, T0
//     * @param saMax   number of iterations before update the temperature
//     */
//    public Tabu(Problem problem, Random random, int tabuListSize) {
//        super(problem, random, "TS");
//
//        // initializing simulated annealing parameters
//        this.tabuListSize = tabuListSize;
//        this.tabuList = new TabuList(this.tabuListSize);
//    }
//	
//    /**
//	 * Find the non-tabu {@link Solution} with the lowest value.<br>
//	 * This method doesn't use any Aspiration Criteria.
//	 */
//	
//    public Solution findBestNeighbor(List<Solution> neighborsSolutions, final List<Solution> solutionsInTabu) {
//		
//		//sort the neighbors
//		Collections.sort(neighborsSolutions, new Comparator<Solution>() {
//			@Override
//			public int compare(Solution a, Solution b) {
//				return ((Double)a.getTotalCost()).compareTo((Double)b.getTotalCost());
//			}
//		});
//
//    	//remove any neighbor that is in tabu list
//		//CollectionUtils.filterInverse(neighborsSolutions, new Predicate<Solution>() {
//		//	@Override
//		//	public boolean evaluate(Solution neighbor) {
//		//		if (solutionsInTabu.contains(neighbor)) {
//		//			System.out.println(count++);
//		//		}
//		//		return solutionsInTabu.contains(neighbor);
//		//	}
//		//});
//				
//		//get the neighbor with lowest value
//		int i = 0;
//		while (solutionsInTabu.contains(neighborsSolutions.get(i))) i++;
//		return neighborsSolutions.get(i);
//	}
//    /**
//     * Executes the Simulated Annealing.
//     *
//     * @param initialSolution the initial (input) solution.
//     * @param timeLimitMillis the time limit (in milliseconds).
//     * @param maxIters        the maximum number of iterations without improvements to execute.
//     * @param output          output PrintStream for logging purposes.
//     * @return the best solution encountered by the SA.
//     */
//    public Solution run(Solution initialSolution, long timeLimitMillis, long maxIters, PrintStream output) {
//        long finalTimeMillis = System.currentTimeMillis() + timeLimitMillis;
//
//        bestSolution = initialSolution;
//        Solution currentSolution = initialSolution.clone();
//        
//        int currentIteration = 0;
//        while (System.currentTimeMillis() < finalTimeMillis) {
//			//List<Solution> candidateNeighbors = currentSolution.getNeighbors();
//			List<Solution> candidateNeighbors = new LinkedList<Solution>();
//			List<Solution> solutionsInTabu = IteratorUtils.toList(tabuList.iterator());
//			
//			for (int i = 0; i < 500; i++) {
//				Solution neighbor = currentSolution.getRandomNeighbor(random);
//				candidateNeighbors.add(neighbor);
//			}
//			
//			Solution bestNeighborFound = findBestNeighbor(candidateNeighbors, solutionsInTabu);
//			if (bestNeighborFound.getTotalCost() < bestSolution.getTotalCost()) {
//				bestSolution = bestNeighborFound;
//				Util.safePrintStatus(output, nIters, bestSolution, currentSolution, "*");
//			}
//			tabuList.add(currentSolution);
//			currentSolution = bestNeighborFound;
//			tabuList.updateSize(currentIteration, bestSolution);
//            nIters++;
//        }
//
//        return bestSolution;
//    }
//
//    /**
//     * Returns the string representation of this heuristic.
//     *
//     * @return the string representation of this heuristic (with parameters values).
//     */
//    public String toString() {
//        return String.format("Simulated Annealing (size=%d)", this.tabuListSize);
//    }
//}
