package dsbp;

import java.io.IOException;
import java.util.Arrays;

public class Test2 {
	public static void main(String[] args) throws IOException {
		/*
		Problem problem = new Problem("/Users/tdhoang/Dropbox/Research/Taxi dispatching/Datasets/in/50x4_1.txt", 100);
		Solution solution = new Solution(problem);
		solution.read("/Users/tdhoang/Documents/workspace/balancing/data/out/50x4_1_out_cvx.txt");
		System.out.println(solution.getBalancingCost());
		System.out.println(solution.getDrivingCost());
		System.out.println(solution.getTotalCost());
		*/
		/*
		double delta = 30;
		double temperature = 10;
		double rate = 1/FastMath.exp(delta / temperature);
		System.out.println(rate);
		System.out.println(-delta/Math.log(rate));
		*/
		int[] a = {1,2,3,4};
		int[] b = {100,200,300,400};
		System.out.println(Arrays.hashCode(a));
		System.out.println(Arrays.hashCode(b));
	}
}
