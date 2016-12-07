package GCNC;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

public class Utilities {

	/**
	 * Given two feature vectors, computes the Pearson's Correlation Coefficient
	 * 
	 * @param x feature vector x
	 * @param y feature vectr y
	 * @return the Pearson's Correlation Coefficient
	 */
	public static double pearsonsCoeff(double[] x, double[] y) {
		PearsonsCorrelation p = new PearsonsCorrelation();
		double corr = p.correlation(x, y);
		
		return corr;
	}
	
	public static double mean(double[] m) {
	    double sum = 0;
	    for (int i = 0; i < m.length; i++) {
	        sum += m[i];
	    }
	    return sum / m.length;
	}
	
	/**
	 * Given a double value and an array, normalized the value in range [0, 1]
	 * using the softmax scaling methos [Theodoridis and Koutroumbas, 2008]
	 * @param w a double value
	 * @param weights the array
	 * @return
	 */
	public static double softmaxScaling(double w, double[] weights) {
		double softmax = 0;
		
		// mean of weights
		Mean m = new Mean();
		double mean = m.evaluate(weights);
		
		//variance of weights
		Variance v = new Variance();
		double var = v.evaluate(weights);
		
		// diff
		double diff = (w - mean) / var;
		
		// normalized value
		softmax = 1 / (1 + Math.exp(-diff));	
		
		return softmax;
	}
	
	public static double zScoreNorm(double a, double[] A) {
		double z = 0;
		
		// mean of A
//		Mean m = new Mean();
		double m = mean(A);// m.evaluate(A);
		
//		//variance of weights
		StandardDeviation s = new StandardDeviation();
		double std = s.evaluate(A);
//		
		z = (a - m) / std;
		
//		System.out.println("z: " + z);
		
		return m;
	}
	
	/**
	 * Given an Arraylist, converts in to an array
	 * @param list the Arraylist
	 * @return the resulting array
	 */
	public static double[] listToArray(ArrayList<Double> list) {
		double[] valuesToArray = new double[list.size()];
		for (int i = 0; i < list.size(); i++) {
			valuesToArray[i] = list.get(i).doubleValue();
		}
		return valuesToArray;
	}
	
	public static int[] listToArrayInt(ArrayList<Integer> list) {
		int[] valuesToArray = new int[list.size()];
		for (int i = 0; i < list.size(); i++) {
			valuesToArray[i] = list.get(i).intValue();
		}
		return valuesToArray;
	}
	
	/**
	 * method that prints graph in node order
	 * @param graph
	 */
	public static void printGraph(UndirectedGraph graph) {
		System.out.println("print graph");
		for (int i=0; i<graph.size(); i++) {
			System.out.println("************************** node "+i+" **************************");
			System.out.println(graph.edgesFrom(i).size() + " " + graph.edgesFrom(i));
			System.out.println(graph.weightsFrom(i).size() + " " + graph.weightsFrom(i));
		}
	} 
	
	/**
	 * method to write a graph to a file
	 * @param graph The undirected graph
	 * @param fileName The filename
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static void graphToFile(UndirectedGraph graph, String fileName) throws FileNotFoundException, IOException { 
		PrintWriter writer = new PrintWriter(fileName, "UTF-8");
		
		for (int i=0; i<graph.size(); i++) {
			int node = i;
			ArrayList<Integer> edges = new ArrayList<Integer>();
			edges.addAll(graph.edgesFrom(i));
			ArrayList<Double> weights = new ArrayList(graph.weightsFrom(i));
			
			for (int j=0; j<edges.size(); j++) {
				String line = node + "	" + edges.get(j) + "	" + weights.get(j);
				writer.println(line);
			}
		}
		
		writer.close();
	}
	
	public static double[][] sumMatrix(double[][] a) {
		double[][] sum_a = new double[a.length][a.length];
		for (int i=0; i<sum_a.length; i++) {
    		for (int j=0; j<sum_a[i].length; j++) {
    			sum_a[i][j] = 0;
    		}
    	}
			
		for (int i=0; i<a.length; i++) {
			int sum = 0;
			for (int j=0; j<a[i].length; j++) {
				sum += a[i][j];
			}
			sum_a[i][i] = sum;
			System.out.println(sum);
		}
		
		return sum_a;
	}
	
	public static double sumArray(double[] a) {
		double sum = 0;
		for (int i=0; i<a.length; i++) {
			sum += a[i];
		}
		return sum;
	}
	
	public static void print2dArray(double[][] array) {
		for (int i=0; i<array.length; i++) {
			System.out.print(i+" -> ");
			for (int j=0; j<array[i].length; j++) {
				System.out.print(array[i][j]+" ");
			}
			System.out.println();
		}
	}

}
