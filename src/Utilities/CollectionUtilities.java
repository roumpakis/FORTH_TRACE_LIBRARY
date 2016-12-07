package Utilities;

import java.util.ArrayList;
import java.util.HashMap;

public class CollectionUtilities {

	public static double getMaxOccurrence(ArrayList<Double> a) {
		int count = 1, tempCount;
		double popular = a.get(0);
		double temp = 0;
		for (int i = 0; i < (a.size() - 1); i++) {
			temp = a.get(i);
			tempCount = 0;
			for (int j = 1; j < a.size(); j++) {
				if (temp == a.get(j))
					tempCount++;
			}
			if (tempCount > count) {
				popular = temp;
				count = tempCount;
			}
		}
		return popular;
	}
	
	public static HashMap<Integer, ArrayList<Double>> transposeHashMap(HashMap<Integer, ArrayList<Double>> data) {
		HashMap<Integer, ArrayList<Double>> dataT = new HashMap<Integer, ArrayList<Double>>();
		//transpose data to get data for each sensor channel
	    for (int i=0; i<data.get(0).size(); i++){
	    	ArrayList<Double> tmp = new ArrayList<Double>();
	    	for (int j=0; j<data.size(); j++){
	    		tmp.add(data.get(j).get(i));	
	    	}
	    	dataT.put(i, tmp);
	    }
	    return dataT;
	}
	
	/**
	 * Method that merges 2 data structures into 1
	 * 
	 * @param featureSet:
	 *            the structure containing all features
	 * @param featureSet2:
	 *            the structure containing Pairwise Correlation features
	 * @return the concatenated data structure
	 */
	public static HashMap<Integer, ArrayList<HashMap<String, Double>>> mergeStructures(
			HashMap<Integer, ArrayList<HashMap<String, Double>>> featureSet,
			ArrayList<ArrayList<HashMap<String, Double>>> featureSet2) {

		HashMap<Integer, ArrayList<HashMap<String, Double>>> featureSet_final = new HashMap<Integer, ArrayList<HashMap<String, Double>>>();

		for (int i = 0; i < featureSet.size(); i++) {
			ArrayList<HashMap<String, Double>> featuresPerChann = featureSet.get(i);
			ArrayList<HashMap<String, Double>> featuresPerChann2 = featureSet2.get(i);
			if (featuresPerChann2 == null)
				continue;

			ArrayList<HashMap<String, Double>> featuresPerChann_final = new ArrayList<HashMap<String, Double>>();

			for (int ii = 0; ii < featuresPerChann.size(); ii++) {
				HashMap<String, Double> h1 = new HashMap<String, Double>();
				HashMap<String, Double> h2 = new HashMap<String, Double>();
				// System.out.println("s:: "+String.format("%03d", ii));
				h1 = featuresPerChann.get(ii);
				for (int j = 0; j < featuresPerChann2.size(); j++) {
					h2 = featuresPerChann2.get(j);
					// System.out.println("h2:"+h2);
					String s = h2.keySet().toString();

					if (s.contains("s" + String.format("%04d", ii))) {
						// System.out.println("sss"+s.substring(1,14));
						String new_s = s.substring(1, 17);
						if (h2.get(new_s) != null) {
							double v = h2.get(new_s);
							HashMap<String, Double> h = new HashMap<String, Double>();
							h.put(new_s.substring(0, 11), v);
							h1.putAll(h);
						}
					}
				}
				featuresPerChann_final.add(h1);

			}
			featureSet_final.put(i, featuresPerChann_final);
		}

		return featureSet_final;
	}
	
	public static double[] listToArray(ArrayList<Double> list) {
		double[] valuesToArray = new double[list.size()];
		for (int i = 0; i < list.size(); i++) {
			valuesToArray[i] = list.get(i);
		}
		return valuesToArray;
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
