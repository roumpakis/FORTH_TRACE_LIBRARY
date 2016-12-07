package GCNC;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
/*****************************************************************************
 * File: UndirectedGraph.java
 * Author: Keith Schwarz (htiek@cs.stanford.edu)
 *
 * A class representing an undirected graph where each edge has an associated 
 * real-valued length.  Internally, the class is represented by an adjacency
 * list where each edges appears twice - once in the forward direction and
 * once in the reverse.  In fact, this implementation was formed by taking
 * a standard adjacency list and then duplicating the logic to ensure each
 * edge appears twice.
 */
import java.util.*; // For HashMap, HashSet
import java.util.Map.Entry;

public final class UndirectedGraph<T> implements Iterable<T> {
    /* A map from nodes in the graph to sets of outgoing edges.  Each
     * set of edges is represented by a map from edges to doubles.
     */
    private Map<T, Set<T>> mGraph = new HashMap<T, Set<T>>();
    
    /* A map from nodes in the graph to sets of outgoing weights.  Each
     * set of weights is represented by a map from edges to doubles.
     */
    private Map<T, ArrayList<Double>> mGraphWeights = new HashMap<T, ArrayList<Double>>();
    
    private Map<Integer, Set<Integer>> mClusters = new HashMap<Integer, Set<Integer>>();

    private Map<Integer, Integer> mNodesClusters = new HashMap<Integer, Integer>();
    
    private double[][] weightMatrix;
    
    private double[][] sumWeightMatrix;
    

	/**
     * Adds a new node to the graph.  If the node already exists, this
     * function is a no-op.
     *
     * @param node The node to add.
     * @return Whether or not the node was added.
     */
    public boolean addNode(T node) {
        /* If the node already exists, don't do anything. */
        if (mGraph.containsKey(node))
            return false;

        /* Otherwise, add the node with an empty set of outgoing edges. */
        mGraph.put(node, new HashSet<T>());
        
        /* Also, add the node with an empty set of outgoing edge weights. */
        mGraphWeights.put(node, new ArrayList<Double>());
        
        mNodesClusters.put((int) node, new Integer(-1));
        
        return true;
    }

    /**
     * Given a node, returns whether that node exists in the graph.
     *
     * @param The node in question.
     * @return Whether that node eixsts in the graph.
     */
    public boolean nodeExists(T node) {
        return mGraph.containsKey(node);
    }

    public void removeNode(T node) {
    	/* Confirm both endpoints exist. */
        if (!mGraph.containsKey(node))
            throw new NoSuchElementException("Node must be in the graph.");
        
        this.mGraph.remove(node);
    }
    
    /**
     * Given two nodes, adds an arc of that length between those nodes.  If 
     * either endpoint does not exist in the graph, throws a 
     * NoSuchElementException.
     *
     * @param one The first node.
     * @param two The second node.
     * @throws NoSuchElementException If either the start or destination nodes
     *                                do not exist.
     */
    public void addEdge(T one, T two, double w) {
        /* Confirm both endpoints exist. */
        if (!mGraph.containsKey(one) || !mGraph.containsKey(two))
            throw new NoSuchElementException("Both nodes must be in the graph.");

        /* Add the edge in both directions. */
        mGraph.get(one).add(two);
        mGraph.get(two).add(one);
    }
    
    public void addEdgeWeight(T one, T two, double w) {
        /* Confirm both endpoints exist. */
        if (!mGraphWeights.containsKey(one) || !mGraphWeights.containsKey(two))
            throw new NoSuchElementException("Both nodes must be in the graph.");

//        if (one!=two) {
        	/* Add the edge weight in both directions. */
        	mGraphWeights.get(one).add(w);
        	mGraphWeights.get(two).add(w);
//        }
    }

    /**
     * Removes the edge between the indicated endpoints from the graph.  If the
     * edge does not exist, this operation is a no-op.  If either endpoint does
     * not exist, this throws a NoSuchElementException.
     *
     * @param one The start node.
     * @param two The destination node.
     * @throws NoSuchElementException If either node is not in the graph.
     */
    public void removeEdge(T one, T two) {
        /* Confirm both endpoints exist. */
        if (!mGraph.containsKey(one) || !mGraph.containsKey(two))
            throw new NoSuchElementException("Both nodes must be in the graph.");

        /* Remove the edges from both adjacency lists. */
        mGraph.get(one).remove(two);
        mGraph.get(two).remove(one);
    }
    
    public void removeEdgeWeight(T one, T two, int w) {
        /* Confirm both endpoints exist. */
        if (!mGraphWeights.containsKey(one) || !mGraphWeights.containsKey(two))
            throw new NoSuchElementException("Both nodes must be in the graph.");

        mGraphWeights.get(one).set((int)w, (double) 0);
    }
    
    public void removeAllEdgesFrom(T one) {
        /* Confirm both endpoints exist. */
        if (!mGraph.containsKey(one))
            throw new NoSuchElementException("Node must be in the graph.");

        /* Remove the edges from both adjacency lists. */
        mGraph.get(one).clear();
    }
    
    public void removeAllEdgeWeightsFrom(T one) {
        /* Confirm both endpoints exist. */
        if (!mGraphWeights.containsKey(one))
            throw new NoSuchElementException("Nodes must be in the graph.");

        mGraphWeights.get(one).clear();
    }

    /**
     * Given two endpoints, returns whether an edge exists between them.  If
     * either endpoint does not exist in the graph, throws a 
     * NoSuchElementException.
     *
     * @param one The first endpoint.
     * @param two The second endpoint.
     * @return Whether an edge exists between the endpoints.
     * @throws NoSuchElementException If the endpoints are not nodes in the 
     *                                graph.
     */
    public boolean edgeExists(T one, T two) {
        /* Confirm both endpoints exist. */
        if (!mGraph.containsKey(one) || !mGraph.containsKey(two))
            throw new NoSuchElementException("Both nodes must be in the graph.");     
        
        /* Graph is symmetric, so we can just check either endpoint. */
        return mGraph.get(one).contains(two);
    }

    /**
     * Given a node in the graph, returns an immutable view of the edges
     * leaving that node.
     *
     * @param node The node whose edges should be queried.
     * @return An immutable view of the edges leaving that node.
     * @throws NoSuchElementException If the node does not exist.
     */
    public Set<T> edgesFrom(T node) {
        /* Check that the node exists. */
        Set<T> arcs = mGraph.get(node);
        if (arcs == null)
            throw new NoSuchElementException("Source node does not exist.");

        return Collections.unmodifiableSet(arcs);
    }
    
    /**
     * Given a node in the graph, returns an immutable view of the edge weights
     * leaving that node.
     *
     * @param node The node whose edges should be queried.
     * @return An immutable view of the edge weights leaving that node.
     * @throws NoSuchElementException If the node does not exist.
     */
    public ArrayList<Double> weightsFrom(T node) {
        /* Check that the node exists. */
        ArrayList<Double> arcs = new ArrayList<Double>(mGraphWeights.get(node));
        if (arcs == null)
            throw new NoSuchElementException("Source node does not exist.");

        return arcs;
    }
    
    /**
     * method to update the edges for a node in the graph
     * @param node
     * @param edges
     */
    public void updateEdgesForNode(T node, Set<T> edges) {
    	/* Check that the node exists. */
        Set<T> arcs = new HashSet<T>(mGraph.get(node));
        if (arcs == null)
            throw new NoSuchElementException("Source node does not exist.");
        
        mGraph.get(node).clear();
//        mGraph.put(node, new HashSet<T>());
        mGraph.get(node).addAll(edges);
    }
    
    /**
     * method to update the edge weights for a node in the graph
     * @param node
     * @param edges
     */
    public void updateWeightsForNode(T node, ArrayList<Double> weights) {
    	/* Check that the node exists. */
        ArrayList<Double> arcs = new ArrayList<Double>(mGraphWeights.get(node));
        if (arcs == null)
            throw new NoSuchElementException("Source node does not exist.");
        
        mGraphWeights.get(node).clear();
        mGraphWeights.get(node).addAll(weights);
    }
    
    /**
     * method to set the cluster value for a specific node
     * @param node
     * @param cluster
     */
    public void setClusterForNode(int cluster, T node) {
    	mNodesClusters.put((int)node, cluster);
    }
    
    /**
     * method to get the cluster value for a specific node
     * @param node
     */
    public int getClusterNode(T node) {
    	return mNodesClusters.get(node);
    }
    
    /**
     * method to add a node to a specific cluster
     * @param node
     * @param cluster
     */
    public void addNodeToCluster(int cluster, int node) {
    	mClusters.get(cluster).add(node);
    }
    
    /**
     * method to get a node from a specific cluster
     * @param cluster
     */
    public Set<Integer> getNodesInCluster(int cluster) {
    	return mClusters.get(cluster);
    }
    
    /**
     * method to set all cluster values for the graph
     * @param clusters
     */
    public void setAllClusters(ArrayList<Integer> clusters) {
    	Set<Integer> uniqueClusters = new HashSet<Integer>(clusters);
    	ArrayList<Integer> un = new ArrayList<Integer> (uniqueClusters);
    	for (int i=0; i<un.size(); i++) {
    		mClusters.put(un.get(i), new HashSet<Integer>());
    	}
    	
    	for (int i=0; i<clusters.size(); i++) {
    		// add node to cluster map
    		mClusters.get(clusters.get(i)).add(i);
    		
    		// add cluster to node map
    		mNodesClusters.put(i, clusters.get(i));
    	}
    }
    
    public Map<Integer, Set<Integer>> getAllClusters() {
    	return mClusters;
    }
    
    /**
     * method that computes the LaplacianEnergy of a graph object
     * @return
     */
    public double LaplacianEnergy() {
    	double LE = 0;
    	/**
    	 * create weight matrix (W) and sum_weight matrix (X)
    	 */
    	setWeightMatrix2();
    	setSumWeightMatrix2();
    	
    	/**
    	 * compute Laplacian Energy
    	 */
    	double sumX = 0;
    	for (int i=0; i<sumWeightMatrix.length; i++) {
    		sumX += sumWeightMatrix[i][i];
    	}
//    	System.out.println("sumX: "+sumX);
    	
    	double sumW = 0;
    	for (int i=0; i<weightMatrix.length; i++) {
    		for (int j=i+1; i<weightMatrix.length; i++) {
    			double w = weightMatrix[i][j];
    			sumW =+ Math.pow(w, 2);
    		}
    	}
    	sumW = 2 * sumW;
    	LE = sumX + sumW;
    	
    	return LE;
    }

    public void W() {
    	double[][] wMatrix = new double[mGraph.size()][mGraph.size()];
    	for (int i=0; i<wMatrix.length; i++) {
    		for (int j=0; j<wMatrix[i].length; j++) {
    			wMatrix[i][j] = 0;
    		}
    	}
    	int c = 0;
    	for (int i=0; i<mGraph.size(); i++) {
    		Set<T> edges = mGraph.get(i);
    		ArrayList<Double> weights = mGraphWeights.get(i);
    		Iterator it = edges.iterator();
    		while (it.hasNext()) {
    			int ed = (int) it.next();
    			System.out.println(i + " -- " + ed);
    		}
    		c++;
    	}
    }
    
    /**
     * Returns whether a given node is contained in the graph.
     *
     * @param The node to test for inclusion.
     * @return Whether that node is contained in the graph.
     */
    public boolean containsNode(T node) {
        return mGraph.containsKey(node);
    }

    /**
     * Returns an iterator that can traverse the nodes in the graph.
     *
     * @return An iterator that traverses the nodes in the graph.
     */
    public Iterator<T> iterator() {
        return mGraph.keySet().iterator();
    }

    /**
     * Returns the number of nodes in the graph.
     *
     * @return The number of nodes in the graph.
     */
    public int size() {
        return mGraph.size();
    }

    /**
     * Returns whether the graph is empty.
     *
     * @return Whether the graph is empty.
     */
    public boolean isEmpty() {
        return mGraph.isEmpty();
    }

    /**
     * Returns a human-readable representation of the graph.
     *
     * @return A human-readable representation of the graph.
     */
    public String toString() {
        return mGraph.toString();
    }
    
    /**
     * Returns a human-readable representation of the graph weights.
     *
     * @return A human-readable representation of the graph weights.
     */
    public String weightstoString() {
        return mGraphWeights.toString();
    }
    
    /**
	 * Getters and Setters
	 */
    
    public Map<T, Set<T>> getmGraph() {
		return mGraph;
	}
    
    public void setmGraph(Map<T, Set<T>> graph) {
//		this.mGraph = graph;
    	for(Entry<T, Set<T>> entry : graph.entrySet()){
		    this.mGraph.put(entry.getKey(), new HashSet<T>(entry.getValue()));
		}
	}

	public Map<T, ArrayList<Double>> getmGraphWeights() {
		return mGraphWeights;
	}
	
	public void setmGraphWeights(Map<T, ArrayList<Double>> weights) {
//		this.mGraphWeights = weights;
		for(Entry<T, ArrayList<Double>> entry : weights.entrySet()){
		    this.mGraphWeights.put(entry.getKey(), new ArrayList<Double>(entry.getValue()));
		}
	}
	
	public double[][] getWeightMatrix() {
		return weightMatrix;
	}

	public void setWeightMatrix(double[][] weightMatrix) {
		this.weightMatrix = weightMatrix;
	}
	
	public void setWeightMatrix2() {
		this.weightMatrix = new double[mGraph.size()][mGraph.size()];
		for (int i=0; i<weightMatrix.length; i++) {
    		for (int j=0; j<weightMatrix[i].length; j++) {
    			weightMatrix[i][j] = 0;
    		}
    	}
		
		for (int i=0; i<weightMatrix.length; i++) {
			Set<Integer> edges = (Set<Integer>) mGraph.get(i);
			ArrayList<Double> weights = mGraphWeights.get(i);
//			System.out.println(i + edges.toString());
//			System.out.println(i + weights.toString());
			Iterator it = edges.iterator();
			int c = 0;
			while (it.hasNext()) { //j
				int e = (int) it.next();
				double w = weights.get(c);
				
				weightMatrix[i][e] = w;
				
				c++;
			}
		}
	}
	
	public void setSumWeightMatrix2() {
		this.sumWeightMatrix = new double[mGraph.size()][mGraph.size()];
		for (int i=0; i<sumWeightMatrix.length; i++) {
    		for (int j=0; j<sumWeightMatrix[i].length; j++) {
    			sumWeightMatrix[i][j] = 0;
    		}
    	}
		
		for (int i=0; i<weightMatrix.length; i++) {
			double sum = 0;
			for (int j=0; j<weightMatrix[i].length; j++) {
//				sum_a[i][j] = 0;
				sum += weightMatrix[i][j];
			}
			sumWeightMatrix[i][i] = sum;
//			System.out.println(sum);
		}
	}

	public double[][] getSumWeightMatrix() {
		return sumWeightMatrix;
	}

	public void setSumWeightMatrix(double[][] sumWeightMatrix) {
		this.sumWeightMatrix = sumWeightMatrix;
	}
	
}