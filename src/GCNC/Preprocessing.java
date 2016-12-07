package GCNC;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.commons.math3.util.MultidimensionalCounter.Iterator;

import Louvain.ModularityOptimizer;
import Utilities.CollectionUtilities;


public class Preprocessing {

	private HashMap<String, ArrayList<Double>> featureSet;
	private String[] featureNames;
	private UndirectedGraph graph;
	private double theta;
	private ArrayList<Integer> clusters;
	private double[][] wMatrix;
	private double[][] xMatrix;

	public Preprocessing(HashMap<String, ArrayList<Double>> featureSet, String[] featureNames, double theta) {
		super();
		this.featureSet = featureSet;
		this.featureNames = featureNames;
		this.graph = new UndirectedGraph();
		this.theta = theta;
		this.clusters = new ArrayList<Integer>();
		
		this.wMatrix = new double[featureSet.size()][featureSet.size()];
		for (int i=0; i<wMatrix.length; i++) {
    		for (int j=0; j<wMatrix[i].length; j++) {
    			wMatrix[i][j] = 0;
    		}
    	}
		
		this.xMatrix = new double[featureSet.size()][featureSet.size()];
	}

	/**
	 * Initializes an UndirectedGraph containing all features as nodes and their 
	 * Pearson's Correlation Coefficient as edges
	 * 
	 * @return
	 */
	public UndirectedGraph initializeGraph(){
		/*
		 * create initial graph with nodes - the features
		 * and edges - their pearsons coefficient
		 */
		for (int i=0; i<featureNames.length; i++){
			// add a node in graph for feature i
			this.graph.addNode(i);
			
			for (int j=0; j<featureNames.length; j++){
				
				// get features i and j
				ArrayList<Double> feat_i = this.featureSet.get(featureNames[i]);
				ArrayList<Double> feat_j = this.featureSet.get(featureNames[j]);
				
				// convert to double[] for correlation computation
				double[] feat_i_arr = Utilities.listToArray(feat_i);
				double[] feat_j_arr = Utilities.listToArray(feat_j);
				
				// compute Pearsons Coefficient between two feature vectors
				double corr = Utilities.pearsonsCoeff(feat_i_arr, feat_j_arr);
				if (i==j) corr = 0;
				
				// add a node in graph for feature j
				this.graph.addNode(j);
				
				// add a graph edge for the two features
				this.graph.addEdge(i, j, corr);
				this.graph.addEdgeWeight(i, j, corr);
			}
		}
		
		/*
		 *  loop to update weights to range [0 1] using softmax scaling
		 */
		for (int i=0; i<featureNames.length; i++) {
			// Arraylist and array with edges corresponding to node i
			ArrayList<Integer> edgesFrom_i = new ArrayList<Integer>(this.graph.edgesFrom(i));
			int[] edgesFrom_i_arr = Utilities.listToArrayInt(edgesFrom_i);
						
			// Arraylist and array with weights corresponding to node i
			ArrayList<Double> weightsFrom_i = new ArrayList<Double>(this.graph.weightsFrom(i));
			double[] weightsFrom_i_arr = Utilities.listToArray(weightsFrom_i);
			
			// update weights for node i
			updateWeights(i, weightsFrom_i_arr);
		}
		
		/*
		 * loop to update edges and weights using the theta threshold
		 * if an edge has weight <= theta it is removed from the graph
		 */
		for (int i=0; i<featureNames.length; i++) {
			// Arraylist and array with edges corresponding to node i
			ArrayList<Integer> edgesFrom_i = new ArrayList<Integer>(this.graph.edgesFrom(i));
			int[] edgesFrom_i_arr = Utilities.listToArrayInt(edgesFrom_i);

			// Arraylist and array with weights corresponding to node i
			ArrayList<Double> weightsFrom_i = new ArrayList<Double>(this.graph.weightsFrom(i));
			double[] weightsFrom_i_arr = Utilities.listToArray(weightsFrom_i);
			
//			double theta = 0.2;
			removeEdgesTheta(i, edgesFrom_i_arr, weightsFrom_i_arr, this.theta);
			
			weightsFrom_i = new ArrayList<Double>(this.graph.weightsFrom(i));
			weightsFrom_i_arr = Utilities.listToArray(weightsFrom_i);
						
		}
		System.out.println("print wMatrix after theta");
		
		return this.graph;
	}

	/**
	 * Given node i and its corresponding weights,
	 * scales each weight value int the range [0,1] using softmax scaling
	 * @param i node 
	 * @param weightsFrom_i_arr corresponding weights
	 */
	public void updateWeights(int i, double[] weightsFrom_i_arr) {
		// init Arraylist for updated weights
		ArrayList<Double> weightsFrom_i_new = new ArrayList<Double>();
		
		// traverse through weights corresponding to node i and perform softmax scaling foreach value
		for (int k=0; k<weightsFrom_i_arr.length; k++) {
			double softmax = 0;
			// softmax scaling for each value
			if (i!=k)
				softmax = Utilities.softmaxScaling(weightsFrom_i_arr[k], weightsFrom_i_arr);
			// add scaled value to new weights Arraylist
			weightsFrom_i_new.add(softmax);
		}
		// update weights for node
		this.graph.updateWeightsForNode(i, weightsFrom_i_new);
	}
	
	/**
	 * Given node i, its corresponding weights and a threshold theta
	 * removes all edges whose weights are below theta
	 * @param i node 
	 * @param weightsFrom_i_arr corresponding weights
	 * @param theta the threshold
	 */
	public void removeEdgesTheta(int i, int[] edgesFrom_i_arr, double[] weightsFrom_i_arr, double theta) {
		// Arraylist with edges and weights to remove
		ArrayList<Integer> badEdges = new ArrayList<Integer>();
		ArrayList<Integer> badWeights = new ArrayList<Integer>();
		
		// Arraylist with edges and weights to keep
		Set<Integer> goodEdges = new HashSet<Integer>();
		ArrayList<Double> goodWeights = new ArrayList<Double>();
		
		// traverse through weights corresponding to node i
		for (int k=0; k<edgesFrom_i_arr.length; k++) {
			int e = edgesFrom_i_arr[k];
			double w = weightsFrom_i_arr[k];
			
			// keep edges with weights greater than theta
			if (w >= theta) {
				goodEdges.add(e);
				goodWeights.add(w);
			}
		} 
		
		// update edges and weights for node i
		this.graph.updateWeightsForNode(i, goodWeights);
		this.graph.updateEdgesForNode(i, goodEdges);
	}
	
	/**
	 * method that performs Louvain Community Detection algorithm on graph
	 * @return the list of clusters for the graph
	 */
	public ArrayList<Integer> louvainClustering() {
		ModularityOptimizer modularity = new ModularityOptimizer();
		
		// init args
		String arguments[] = new String[] {"files/features_graph.txt", 			/* inputFileName */
											"files/features_communities.txt", 	/* outputFileName */
											"1", 								/* modularityFunction */
											"1.0", 								/* resolution */
											"1", 								/* algorithm */
											"10", 								/* nRandomStarts */
											"10", 								/* nIterations */
											"0", 								/* randomSeed */
											"0" 								/* printOutput */
		};

		// call louvain clustering
		try {
			this.clusters = ModularityOptimizer.runOptimizer(arguments);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// set graph clusters
		this.graph.setAllClusters(clusters);
		
		return this.clusters;
	}

	
	/**
	 * method that writes graph to file in order to use it as an input in ModularityOptimizer
	 * @param graph
	 * @param fileName
	 */
	public void writeGraphToFile(UndirectedGraph graph, String fileName) {
		try {
			Utilities.graphToFile(graph, fileName);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * Getters and Setters
	 */
	
	public UndirectedGraph getGraph() {
		return graph;
	}

	public void setGraph(UndirectedGraph graph) {
		this.graph = graph;
	}
	
	public double getTheta() {
		return theta;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}
}
