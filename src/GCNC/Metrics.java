package GCNC;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

public class Metrics {

	/**
	 * Given an Arraylist, this method computes its term variance
	 * @param A The Arraylist
	 * @return The TV value
	 */
	public static double TV(ArrayList<Double> A) {
		double tv = 0;
		double[] A_arr = Utilities.listToArray(A);
		int nPatterns = A.size();
		
		double sum_diff = 0;
		for (int i=0; i<nPatterns; i++) {
			double a = A.get(i);
			// mean of feature vector
			double mean_A = Utilities.mean(A_arr);
			
			// diff
			double diff = Math.pow((a - mean_A), 2);

			// compute the sum
			sum_diff += diff;
		}
		
		// compute term variance value
		tv = sum_diff / nPatterns;
		
		return tv;
	}
	
	/**
	 * Method that computes the Laplacian Centrality of a node in its own cluster
	 * @param <T>
	 * @param node
	 * @param nodeValues
	 * @param cluster
	 * @return
	 */
	public static <T> double LC(int node, ArrayList<Double> nodeValues, Set<Integer> cluster, UndirectedGraph graph) {	
		double LC = 0;
		UndirectedGraph tmpGraph = new UndirectedGraph();
		
		Map<T, Set<T>> mGraph = new HashMap<T, Set<T>>();
		mGraph = graph.getmGraph();
		Map<T, ArrayList<Double>> mGraphWeights = new HashMap<T, ArrayList<Double>>();
		mGraphWeights = graph.getmGraphWeights();
		
		tmpGraph.setmGraph(mGraph);
		tmpGraph.setmGraphWeights(mGraphWeights);
		
		
		Iterator it = cluster.iterator();
		while (it.hasNext()) {
			int i = (int) it.next();
			tmpGraph.updateEdgesForNode(i, graph.edgesFrom(i));
			tmpGraph.updateWeightsForNode(i, graph.weightsFrom(i));
		}
		
		double LE = 0;
		LE = tmpGraph.LaplacianEnergy();
		
		/**
		 * compute LE of the new graph after removing the i-th feature
		 */
		tmpGraph.removeAllEdgesFrom(node);
		tmpGraph.removeAllEdgeWeightsFrom(node);
		double LEtmp = 0;
		LEtmp = tmpGraph.LaplacianEnergy();
		
		/**
		 * compute LC for cluster
		 */
		if (LE == 0)
			LC = 0;
		else
			LC = (LE - LEtmp)/LE;

		return LC;
	}
}
