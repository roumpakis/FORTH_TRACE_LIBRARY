import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeSet;

import Acquisition.AcquisitionMultiLocation;
import Acquisition.AcquisitionSingleLocation;
import FeatureExtraction.FeatureExtraction;
import FeatureSelection.FSSA;
import FeatureSelection.ReliefFAttributeEval;
import FeatureSelection.SequentialFloatingFS;
import GCNC.GCNC;
import GCNC.Preprocessing;
import GCNC.UndirectedGraph;
import Louvain.ModularityOptimizer;
import Segmentation.MultiLocationData;
import Segmentation.SegmentationMultiLocation;
import Segmentation.SegmentationSingleLocation;
import Utilities.*;
import net.sf.javaml.utils.ArrayUtils;
import weka.core.Instances;

public class Main {
	public static void main(String[] args) {
		int mode = 1;
		
		HAPT_experiment();
		TRACE_experiment();
		TRACE_experiment_per_device();
	}

	public static void HAPT_experiment(){
		String path = "/path/to/HAPT/";
		
		/** 
		 * initialize 10 random participants
		 */
		ArrayList<Integer> participants = new ArrayList<Integer>();
		int min = 1, max = 14;
		Random rn = new Random();
		for (int i=0; i<10; i++){
			int participant_no = rn.nextInt(max - min + 1) + min;
			while (participants.contains(participant_no))
				participant_no = rn.nextInt(max - min + 1) + min;
			participants.add(participant_no);
		}
		
		ArrayList<Results> results = new ArrayList<Results>();
		for (int p=0; p<10; p++){ // start 10 random participants
			int participant_no = participants.get(p);
			String directory = "user"+String.format("%02d", participant_no);
			String trial_no = "exp01";
			String filename = path+trial_no+"_"+directory+".csv";
			System.out.println(filename);
			
			// **** DATA ACQUISITION **** //
			AcquisitionSingleLocation acquisition = new AcquisitionSingleLocation(filename);
			acquisition.parseDataSingleLocation();
			HashMap<Integer, ArrayList<Double>> rawData = acquisition.getDataT();
			ArrayList<Double> labels = acquisition.getLabels();
			
			// **** SEGMENTATION **** //
			int[] windowSizeSec = {2, 5, 10, 20};
			double overlapping_a = 0.5;
			
			for (int w=0; w<windowSizeSec.length; w++) {
				SegmentationSingleLocation segmentation = new SegmentationSingleLocation(rawData, labels);
				segmentation.SegmentationOverlapping(windowSizeSec[w], overlapping_a);
				HashMap<Integer, ArrayList<ArrayList<Double>>> segmentedData = segmentation.getSegmentedData();
				ArrayList<Double> finalLabels = segmentation.getFinalLabels();
	
				// **** FEATURE EXTRACTION **** //
				boolean normalized = true;
				long start_time_FE = System.nanoTime();
				FeatureExtraction FE = new FeatureExtraction(segmentedData, finalLabels, 0, windowSizeSec[w]);
				FE.createFeatureSet(normalized);
				long elapsed_time_FE = System.nanoTime() - start_time_FE;
				
				HashMap<Integer, ArrayList<HashMap<String,Double>>> featureSet = FE.getInitialFeatureSet();
				HashMap<String, ArrayList<Double>> finalFeatureSet_all = FE.getFinalFeatureSet();
				HashMap<String, ArrayList<Double>> finalFeatureSet = new HashMap<String, ArrayList<Double>>();
				finalFeatureSet.putAll(finalFeatureSet_all);
				
				if (normalized) {
					ArrayList<String> badFeatures = FE.getBadFeatures();
					System.out.println("badFeatures: "+badFeatures);
					for (int i=0; i<badFeatures.size(); i++) {
						finalFeatureSet.remove(badFeatures.get(i));
					}
				}
				System.out.println("prev size: "+finalFeatureSet_all.size()+" curr size: "+finalFeatureSet.size());
				
				TreeSet<String> featureNames = new TreeSet<String>(finalFeatureSet.keySet());
				String[] featureNames_Arr = featureNames.toArray(new String[featureNames.size()]);
				
				TreeSet<String> featureNames_all = new TreeSet<String>(finalFeatureSet_all.keySet());
				String[] featureNames_Arr_all = featureNames_all.toArray(new String[featureNames_all.size()]);
				
				
				// **** FEATURE SELECTION **** //
				/**
				 * SFFS
				 */
				SequentialFloatingFS sffs = new SequentialFloatingFS(finalFeatureSet, finalLabels);
				
				long start_time_sffs = System.nanoTime();
				HashMap<String, ArrayList<Double>> sffs_subset = sffs.SFFS(false);
				long elapsed_time_sffs = System.nanoTime() - start_time_sffs;
				System.out.println("elapsed time SFFS: "+elapsed_time_sffs);
				
				System.out.println(sffs_subset.keySet());

				ArrayList<String> SFFSdominantFeaturesStr = new ArrayList<String>(sffs_subset.keySet());
				ArrayList<Integer> SFFSdominantFeaturesInd = new ArrayList<Integer>();
				int SFFS_subset_size = SFFSdominantFeaturesStr.size();
				
				for (int i=0; i<SFFSdominantFeaturesStr.size(); i++){ 
					for (int j=0; j<featureNames_Arr.length; j++){
						if (SFFSdominantFeaturesStr.get(i).equals(featureNames_Arr[j]))
							SFFSdominantFeaturesInd.add(j);
					}
				}
				
				Results res0 = new Results(participant_no, windowSizeSec[w]);
				res0.setAlgNo("0");
				res0.setWinsizeSec(windowSizeSec[w]);
				res0.setFeaturesNo(sffs_subset.size());
				if (res0.getFeaturesNo() > 0)
					res0.setDominantFeatures(SFFSdominantFeaturesInd);
				res0.setExecTimeFS(elapsed_time_sffs);
				res0.setExecTimeFE(elapsed_time_FE);
				
				/**
				 * FSSA
				 */
				FSSA fssa = new FSSA(finalFeatureSet, featureNames_Arr);
				ArrayList<Results> results1 = new ArrayList<Results>();
				double[] percentages = new double[] {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
				for (int perc=0; perc<percentages.length; perc++) {
					double k =  percentages[perc] * finalFeatureSet.size();
					System.out.println("k"+percentages[perc]);
					long start_time_fssa = System.nanoTime();
					HashMap<String, Double> fssa_subset = fssa.FSSA((int) k);
					long elapsed_time_fssa = System.nanoTime() - start_time_fssa;
					System.out.println("elapsed_time_fssa: "+elapsed_time_fssa);
					
					System.out.println("FSSA subset: "+fssa_subset.keySet());
					ArrayList<String> FSSAdominantFeaturesStr = new ArrayList<String>(fssa_subset.keySet());
					ArrayList<Integer> FSSAdominantFeaturesInd = new ArrayList<Integer>();
					int fssa_subset_size = fssa_subset.keySet().size();
					
					for (int i=0; i<FSSAdominantFeaturesStr.size(); i++){ 
						for (int j=0; j<featureNames_Arr_all.length; j++){
							if (FSSAdominantFeaturesStr.get(i).equals(featureNames_Arr_all[j]))
								FSSAdominantFeaturesInd.add(j);
						}
					}
					
					Results res1 = new Results(participant_no, windowSizeSec[w]);
					res1.setAlgNo("1-"+percentages[perc]*100);
					res1.setWinsizeSec(windowSizeSec[w]);
					res1.setFeaturesNo(fssa_subset_size);
					if (res1.getFeaturesNo() > 0)
						res1.setDominantFeatures(FSSAdominantFeaturesInd);
					res1.setExecTimeFS(elapsed_time_fssa);
					res1.setExecTimeFE(elapsed_time_FE);
					
					results1.add(res1);
				}
				
				
				/**
				 * RELIEF
				 */
				ReliefFAttributeEval relief = new ReliefFAttributeEval();
				Instances ReliefInput = relief.createReliefInput(fssa.getFSSA_input(), featureNames_Arr, finalLabels);
				ArrayList<Integer> ReliefSubset = new ArrayList<Integer>();
				Results res2 = new Results(participant_no, windowSizeSec[w]);
				
				try {
					String[] options = {"-M", "-1", "-D", "1", "-K", "10", "-W", "-A", "2"};
					relief.setOptions(options);
					relief.setSampleSize(-1);
					relief.setSeed(1);
					relief.setNumNeighbours(10);
					long start_time_relief = System.nanoTime();
					relief.buildEvaluator(ReliefInput);
					long elapsed_time_relief = System.nanoTime() - start_time_relief;
					System.out.println("elapsed_time_relief: "+elapsed_time_relief);
					
					
					ReliefSubset = relief.getIndexes();
					System.out.println("ReliefSubset"+ReliefSubset);
					int Relief_subset_size = ReliefSubset.size();
					res2.setAlgNo("2");
					res2.setWinsizeSec(windowSizeSec[w]);
					res2.setFeaturesNo(Relief_subset_size);
					
					if (!relief.containsNan())
						res2.setDominantFeatures(ReliefSubset);
					else {
						res2.setDominantFeatures(new ArrayList<Integer>());
						res2.setFeaturesNo(0);
					}
					res2.setExecTimeFS(elapsed_time_relief);
					res2.setExecTimeFE(elapsed_time_FE);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				/**
				 * GCNC
				 */
				System.out.println("GCNC part"+participant_no);
				int nPatterns = finalLabels.size();
				System.out.println("nPatterns: "+nPatterns);	
				
				/**
				 * Initialize GCNC
				 * PREPROCESSING
				 */
				double theta = 0.5;
				Preprocessing prep = new Preprocessing(finalFeatureSet, featureNames_Arr, theta);
				UndirectedGraph graph =  prep.initializeGraph();
				prep.writeGraphToFile(graph, "files/features_graph.txt");	
				ArrayList<Integer> graphClusters = prep.louvainClustering();
				System.out.println("graphClusters: "+graphClusters);
				System.out.println(graph.getAllClusters().entrySet());
				
				/**
				 * Start GCNC
				 */
				ArrayList<Results> results4 = new ArrayList<Results>();
				
				double delta = 0.0;
				double om = 0;
				GCNC gcnc= new GCNC(finalFeatureSet, featureNames_Arr, nPatterns, graph, delta, om);
				long start_time_GCNC = System.nanoTime();
				gcnc.core();
				long elapsed_time_GCNC = System.nanoTime() - start_time_GCNC;
				System.out.println("elapsed_time_GCNC: "+elapsed_time_GCNC);
				
				ArrayList<Integer> GCNC_dominantFeatInd = gcnc.getReducedSet();
				
				Results res3 = new Results(participant_no, windowSizeSec[w]);
				res3.setAlgNo("3");
				res3.setWinsizeSec(windowSizeSec[w]);
				res3.setFeaturesNo(GCNC_dominantFeatInd.size());
				if (res3.getFeaturesNo() > 0)
					res3.setDominantFeatures(GCNC_dominantFeatInd);
				res3.setExecTimeFS(elapsed_time_GCNC);
				res3.setExecTimeFE(elapsed_time_FE);
			
				results4.add(res3);
				
				
				results.add(res0);
				results.addAll(results1);
				results.add(res2);
				
				results.addAll(results4);
				Utilities.IO.writeFeatureSetToFile("files/HAPT_experiment/part"+participant_no+"winsize"+windowSizeSec[w]+"features.csv", finalFeatureSet_all,
											finalLabels, featureNames_Arr_all);
			}
			Utilities.IO.writeResultsToFile("files/HAPT_experiment/stats.csv", results);
		}
	}
	
	/**
	 * experiment on TRACE dataset
	 * concatenated files for each participant
	 */
	public static void TRACE_experiment(){
		String path = "/path/to/TRACE/";
		
		/** 
		 * initialize 10 random participants
		 */
		ArrayList<Integer> participants = new ArrayList<Integer>();
		int [] parts = new int[] {0, 1, 2, 3, 5, 6, 7, 9, 10, 12, 13, 14};
		int min = 0, max = 11;
		Random rn = new Random();
		for (int i=0; i<10; i++){
			int participant_no = parts[rn.nextInt(max - min + 1) + min];
			while ( (participants.contains(participant_no)) )
				participant_no = parts[rn.nextInt(max - min + 1) + min];
			participants.add(participant_no);
		}
		System.out.println(participants);
		
		for (int p=0; p<10; p++){ // start 10 random participants
			ArrayList<Results> results = new ArrayList<Results>();
			int participant_no = participants.get(p);
			String filename = path+"alldata_concatenated_part"+participant_no+".csv";
			System.out.println(filename);
			
			// **** DATA ACQUISITION **** //
			AcquisitionMultiLocation acquisition = new AcquisitionMultiLocation(filename);
			acquisition.parseDataMultiLocation();
			HashMap<Integer, ArrayList<Double>> rawData = acquisition.getData();
			HashMap<Integer, ArrayList<Double>> rawDataT = acquisition.getDataT();
			ArrayList<Double> labels = acquisition.getLabels();
			ArrayList<Double> timestamps = acquisition.getTimestamps();
			ArrayList<Double> device_ids = acquisition.getDevice_ids();
			
			// **** SEGMENTATION **** //
			double overlapping_a = 0.5;
			int[] windowSizeSec = {2,5,10,20};
			
			for (int w=0; w<windowSizeSec.length; w++){ // start window size experiments
				System.out.println("participant: "+ participant_no+" winsize: "+windowSizeSec[w]);
				SegmentationMultiLocation segmentation = new SegmentationMultiLocation(rawDataT, labels, timestamps, device_ids);
				segmentation.SegmentationOverlapping(windowSizeSec[w], overlapping_a);
				HashMap<Integer, ArrayList<ArrayList<Double>>> segmentedData = segmentation.getSegmentedData();
				ArrayList<ArrayList<Double>> segmentedDevice_ids = segmentation.getSegmentedDevice_ids();
				ArrayList<Double> finalLabels = segmentation.getFinalLabels();
				System.out.println("no segments:"+finalLabels.size());
				
				MultiLocationData all_Data = new MultiLocationData(segmentedData, finalLabels, segmentedDevice_ids);
				all_Data.separateData();
				HashMap<Integer, HashMap<Integer, ArrayList<ArrayList<Double>>>> separated_data = all_Data.getSeparated_data();
				ArrayList<ArrayList<Double>> separated_device_ids = all_Data.getAll_device_ids();
				ArrayList<Integer> badIndexes = all_Data.getBadIndexes();
				System.out.println("bad"+badIndexes);
				for (int ind=0; ind<badIndexes.size(); ind++) {
					finalLabels.remove(badIndexes.get(ind).intValue());
				}
				
				// **** FEATURE EXTRACTION **** //
				HashMap<String, ArrayList<Double>> featureSetAllDevices = new HashMap<String, ArrayList<Double>>();
				HashMap<String, ArrayList<Double>> featureSetAllDevices2 = new HashMap<String, ArrayList<Double>>();
				ArrayList<Integer> badSegments = new ArrayList<Integer>();
				long start_time_FE = System.nanoTime();
				for (int dev=1; dev<=5; dev++) {
					System.out.println("Feature Extraction dev" +dev);
					// ** extract features for each separate device ** //
					boolean normalized = true;
					FeatureExtraction FE = new FeatureExtraction(
							separated_data.get(dev), finalLabels, 
							separated_device_ids, dev, windowSizeSec[w]);
					FE.createFeatureSet(normalized);
					HashMap<Integer, ArrayList<HashMap<String,Double>>> featureSet = FE.getInitialFeatureSet();
					HashMap<String, ArrayList<Double>> finalFeatureSet_all = FE.getFinalFeatureSet();
					badSegments = FE.getBadSegments();
					featureSetAllDevices2.putAll(finalFeatureSet_all);
					
					HashMap<String, ArrayList<Double>> finalFeatureSet = new HashMap<String, ArrayList<Double>>();
					finalFeatureSet.putAll(finalFeatureSet_all);
					
					// remove bad features
					if (normalized) {
						ArrayList<String> badFeatures = FE.getBadFeatures();
						for (int i=badFeatures.size()-1; i>=0; i--) {
							finalFeatureSet.remove(badFeatures.get(i));
						}
					}
					featureSetAllDevices.putAll(finalFeatureSet);
				}
				long elapsed_time_FE = System.nanoTime() - start_time_FE;
				System.out.println("elapsed_time_FE: "+elapsed_time_FE);
				
				//remove bad segments from labels as well
				for (int l=0; l<badSegments.size(); l++) {
					finalLabels.set(badSegments.get(l), Double.NaN);
				}
				finalLabels.removeAll(Collections.singleton(Double.NaN));
				System.out.println("no segments final:"+finalLabels.size());
				TreeSet<String> featureNames = new TreeSet<String>(featureSetAllDevices.keySet());
				String[] featureNames_Arr = featureNames.toArray(new String[featureNames.size()]);
				
				TreeSet<String> featureNames2 = new TreeSet<String>(featureSetAllDevices2.keySet());
				String[] featureNames_Arr2 = featureNames2.toArray(new String[featureNames2.size()]);
				
				// **** FEATURE SELECTION **** //
				/**
				 * SFFS
				 */
				SequentialFloatingFS sffs = new SequentialFloatingFS(featureSetAllDevices, finalLabels);
				long start_time_sffs = System.nanoTime();
				HashMap<String, ArrayList<Double>> sffs_subset = sffs.SFFS(false);
				long elapsed_time_sffs = System.nanoTime() - start_time_sffs;
				System.out.println("elapsed time SFFS: "+elapsed_time_sffs);
				System.out.println("SFFS subset: "+sffs_subset.keySet());
				
				TreeSet<String> featureNames_all = new TreeSet<String>(featureSetAllDevices2.keySet());
				String[] featureNames_Arr_all = featureNames_all.toArray(new String[featureNames_all.size()]);
				
				ArrayList<String> SFFSdominantFeaturesStr = new ArrayList<String>(sffs_subset.keySet());
				ArrayList<Integer> SFFSdominantFeaturesInd = new ArrayList<Integer>();
				int SFFS_subset_size = SFFSdominantFeaturesStr.size();
				
				for (int i=0; i<SFFSdominantFeaturesStr.size(); i++){ 
					for (int j=0; j<featureNames_Arr2.length; j++){
						if (SFFSdominantFeaturesStr.get(i).equals(featureNames_Arr2[j]))
							SFFSdominantFeaturesInd.add(j);
					}
				}
				
				Results res0 = new Results(participant_no, windowSizeSec[w]);
				res0.setAlgNo("0");
				res0.setWinsizeSec(windowSizeSec[w]);
				res0.setFeaturesNo(sffs_subset.size());
				if (res0.getFeaturesNo() > 0)
					res0.setDominantFeatures(SFFSdominantFeaturesInd);
				res0.setExecTimeFS(elapsed_time_sffs);
				res0.setExecTimeFE(elapsed_time_FE);
				
				/**
				 * FSSA
				 */
				FSSA fssa = new FSSA(featureSetAllDevices, featureNames_Arr);
				ArrayList<Results> results1 = new ArrayList<Results>();
				double[] percentages = new double[] {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
				for (int perc=0; perc<percentages.length; perc++) {
					double k =  percentages[perc] * featureSetAllDevices.size();
					System.out.println("k"+percentages[perc]);
					long start_time_fssa = System.nanoTime();
					HashMap<String, Double> fssa_subset = fssa.FSSA((int) k);
					long elapsed_time_fssa = System.nanoTime() - start_time_fssa;
					System.out.println("elapsed_time_fssa: "+elapsed_time_fssa);
					
					System.out.println("FSSA subset: "+fssa_subset.keySet());
					ArrayList<String> FSSAdominantFeaturesStr = new ArrayList<String>(fssa_subset.keySet());
					ArrayList<Integer> FSSAdominantFeaturesInd = new ArrayList<Integer>();
					int fssa_subset_size = fssa_subset.keySet().size();
					
					for (int i=0; i<FSSAdominantFeaturesStr.size(); i++){ 
						for (int j=0; j<featureNames_Arr2.length; j++){
							if (FSSAdominantFeaturesStr.get(i).equals(featureNames_Arr2[j]))
								FSSAdominantFeaturesInd.add(j);
						}
					}
					
					Results res1 = new Results(participant_no, windowSizeSec[w]);
					res1.setAlgNo("1-"+percentages[perc]*100);
					res1.setWinsizeSec(windowSizeSec[w]);
					res1.setFeaturesNo(fssa_subset_size);
					if (res1.getFeaturesNo() > 0)
						res1.setDominantFeatures(FSSAdominantFeaturesInd);
					res1.setExecTimeFS(elapsed_time_fssa);
					res1.setExecTimeFE(elapsed_time_FE);
					
					results1.add(res1);
				}
		
				
				/**
				 * RELIEF
				 */
				ReliefFAttributeEval relief = new ReliefFAttributeEval();
				Instances ReliefInput = relief.createReliefInput(fssa.getFSSA_input(), featureNames_Arr, finalLabels);
				ArrayList<Integer> ReliefSubset = new ArrayList<Integer>();
				ArrayList<String> ReliefSubsetStr = new ArrayList<String>();
				ArrayList<Integer> ReliefSubsetInd = new ArrayList<Integer>();
				Results res2 = new Results(participant_no, windowSizeSec[w]);
				
				try {
					String[] options = {"-M", "-1", "-D", "1", "-K", "10", "-W", "-A", "2"};
					relief.setOptions(options);
					relief.setSampleSize(-1);
					relief.setSeed(1);
					relief.setNumNeighbours(10);
					long start_time_relief = System.nanoTime();
					relief.buildEvaluator(ReliefInput);
					long elapsed_time_relief = System.nanoTime() - start_time_relief;
					System.out.println("elapsed_time_relief: "+elapsed_time_relief);
					
					ReliefSubset = relief.getIndexes();
//					System.out.println("ReliefSubset"+ReliefSubset);
					ReliefSubsetStr = relief.getReliefSubset();
					System.out.println("ReliefSubsetStr"+ReliefSubsetStr);
					
					for (int i=0; i<ReliefSubsetStr.size(); i++){ 
						for (int j=0; j<featureNames_Arr2.length; j++){
							if (ReliefSubsetStr.get(i).equals(featureNames_Arr2[j]))
								ReliefSubsetInd.add(j);
						}
					}
					System.out.println("ReliefSubsetInd"+ReliefSubsetInd);
					int Relief_subset_size = ReliefSubset.size();
					res2.setAlgNo("2");
					res2.setWinsizeSec(windowSizeSec[w]);
					res2.setFeaturesNo(Relief_subset_size);
					
					if (!relief.containsNan())
						res2.setDominantFeatures(ReliefSubsetInd);
					else {
						res2.setDominantFeatures(new ArrayList<Integer>());
						res2.setFeaturesNo(0);
					}
					res2.setExecTimeFS(elapsed_time_relief);
					res2.setExecTimeFE(elapsed_time_FE);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				/**
				 * GCNC
				 */
				System.out.println("GCNC part"+participant_no);
				int nPatterns = finalLabels.size();
				System.out.println("nPatterns: "+nPatterns);	
				
				
				double theta = 0.5;
				Preprocessing prep = new Preprocessing(featureSetAllDevices, featureNames_Arr, theta);
				UndirectedGraph graph =  prep.initializeGraph();
//				p.printGraph();
				prep.writeGraphToFile(graph, "files/features_graph.txt");	
				ArrayList<Integer> graphClusters = prep.louvainClustering();
				System.out.println("graphClusters: "+graphClusters);
				System.out.println(graph.getAllClusters().entrySet());
				
				
				ArrayList<Results> results4 = new ArrayList<Results>();
				
				double delta = 0.0;
				double om = 0;// Math.floor(featureNames_Arr.length/2);
				GCNC gcnc= new GCNC(featureSetAllDevices, featureNames_Arr, nPatterns, graph, delta, om);
				long start_time_GCNC = System.nanoTime();
				gcnc.core();
				long elapsed_time_GCNC = System.nanoTime() - start_time_GCNC;
				System.out.println("elapsed_time_GCNC: "+elapsed_time_GCNC);
				
				ArrayList<Integer> GCNC_dominantFeatInd = gcnc.getReducedSet();
				
				Results res3 = new Results(participant_no, windowSizeSec[w]);
				res3.setAlgNo("3");
				res3.setWinsizeSec(windowSizeSec[w]);
				res3.setFeaturesNo(GCNC_dominantFeatInd.size());
				if (res3.getFeaturesNo() > 0)
					res3.setDominantFeatures(GCNC_dominantFeatInd);
				res3.setExecTimeFS(elapsed_time_GCNC);
				res3.setExecTimeFE(elapsed_time_FE);				
				results4.add(res3);
				
				
				
				results.add(res0);
				results.addAll(results1);
				results.add(res2);

				results.addAll(results4);
				
				Utilities.IO.writeFeatureSetToFile("files/TRACE_experiment/part"+participant_no+"winsize"+windowSizeSec[w]+"features.csv", featureSetAllDevices2,
											finalLabels, featureNames_Arr_all);
			}
			Utilities.IO.writeResultsToFile("files/TRACE_experiment/stats.csv", results);
		}
		
	}

	/**
	 * experiment on TRACE dataset
	 * for each device id separately
	 */
	public static void TRACE_experiment_per_device(){
		String path = "/path/to/TRACE/";
		
		/** 
		 * initialize 10 random participants
		 */
		ArrayList<Integer> participants = new ArrayList<Integer>();
		int min = 1, max = 14;
		Random rn = new Random();
		for (int i=0; i<10; i++){
			int participant_no = rn.nextInt(max - min + 1) + min;
			while (participants.contains(participant_no))
				participant_no = rn.nextInt(max - min + 1) + min;
			participants.add(participant_no);
		}
		
		/**
		 * 5 devices
		 */
		for (int dev=1; dev<=5; dev++) {
			ArrayList<Results> results = new ArrayList<Results>();
			
			/**
			 * 10 participants
			 */
			for (int p=0; p<10; p++){ // start 10 random participants

				int participant_no = participants.get(p);
				String filename = path+"part"+participant_no+"/part"+participant_no+"dev"+dev+".csv";

				System.out.println(filename);
				
				// **** DATA ACQUISITION **** //
				AcquisitionSingleLocation acquisition = new AcquisitionSingleLocation(filename);
				acquisition.parseDataSingleLocation();
				HashMap<Integer, ArrayList<Double>> rawData = acquisition.getDataT();
				ArrayList<Double> labels = acquisition.getLabels();
				
				// **** SEGMENTATION **** //
				int[] windowSizeSec = {2, 5, 10, 20};
				double overlapping_a = 0.5;
				
				for (int w=0; w<windowSizeSec.length; w++) {
					System.out.println("winsize: "+windowSizeSec[w]);
					
					SegmentationSingleLocation segmentation = new SegmentationSingleLocation(rawData, labels);
					segmentation.SegmentationOverlapping(windowSizeSec[w], overlapping_a);
					HashMap<Integer, ArrayList<ArrayList<Double>>> segmentedData = segmentation.getSegmentedData();
					ArrayList<Double> finalLabels = segmentation.getFinalLabels();
		
					// **** FEATURE EXTRACTION **** //
					boolean normalized = true;
					long start_time_FE = System.nanoTime();
					FeatureExtraction FE = new FeatureExtraction(segmentedData, finalLabels, 0, windowSizeSec[w]);
					FE.createFeatureSet(normalized);
					long elapsed_time_FE = System.nanoTime() - start_time_FE;
					
					HashMap<Integer, ArrayList<HashMap<String,Double>>> featureSet = FE.getInitialFeatureSet();
					HashMap<String, ArrayList<Double>> finalFeatureSet_all = FE.getFinalFeatureSet();
					HashMap<String, ArrayList<Double>> finalFeatureSet = new HashMap<String, ArrayList<Double>>();
					finalFeatureSet.putAll(finalFeatureSet_all);
					
					if (normalized) {
						ArrayList<String> badFeatures = FE.getBadFeatures();
						System.out.println("badFeatures: "+badFeatures);
						for (int i=0; i<badFeatures.size(); i++) {
							finalFeatureSet.remove(badFeatures.get(i));
						}
					}
					System.out.println("prev size: "+finalFeatureSet_all.size()+" curr size: "+finalFeatureSet.size());
					
					TreeSet<String> featureNames = new TreeSet<String>(finalFeatureSet.keySet());
					String[] featureNames_Arr = featureNames.toArray(new String[featureNames.size()]);
					
					TreeSet<String> featureNames_all = new TreeSet<String>(finalFeatureSet_all.keySet());
					String[] featureNames_Arr_all = featureNames_all.toArray(new String[featureNames_all.size()]);
					
					// **** FEATURE SELECTION **** //
					/**
					 * SFFS
					 */
					SequentialFloatingFS sffs = new SequentialFloatingFS(finalFeatureSet, finalLabels);
					
					long start_time_sffs = System.nanoTime();
					HashMap<String, ArrayList<Double>> sffs_subset = sffs.SFFS(false);
					long elapsed_time_sffs = System.nanoTime() - start_time_sffs;
					System.out.println("elapsed time SFFS: "+elapsed_time_sffs);				    
					System.out.println(sffs_subset.keySet());
					
					
					ArrayList<String> SFFSdominantFeaturesStr = new ArrayList<String>(sffs_subset.keySet());
					ArrayList<Integer> SFFSdominantFeaturesInd = new ArrayList<Integer>();
					int SFFS_subset_size = SFFSdominantFeaturesStr.size();
					
					for (int i=0; i<SFFSdominantFeaturesStr.size(); i++){ 
						for (int j=0; j<featureNames_Arr.length; j++){
							if (SFFSdominantFeaturesStr.get(i).equals(featureNames_Arr[j]))
								SFFSdominantFeaturesInd.add(j);
						}
					}
					
					Results res0 = new Results(participant_no, windowSizeSec[w]);
					res0.setAlgNo("0");
					res0.setWinsizeSec(windowSizeSec[w]);
					res0.setFeaturesNo(sffs_subset.size());
					if (res0.getFeaturesNo() > 0)
						res0.setDominantFeatures(SFFSdominantFeaturesInd);
					res0.setExecTimeFS(elapsed_time_sffs);
					res0.setExecTimeFE(elapsed_time_FE);
					
					/**
					 * FSSA
					 */
					FSSA fssa = new FSSA(finalFeatureSet, featureNames_Arr);
					ArrayList<Results> results1 = new ArrayList<Results>();
					double[] percentages_fssa = new double[] {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
					for (int perc=0; perc<percentages_fssa.length; perc++) {
						double k =  percentages_fssa[perc] * finalFeatureSet.size();
						System.out.println("k"+percentages_fssa[perc]);
						long start_time_fssa = System.nanoTime();
						HashMap<String, Double> fssa_subset = fssa.FSSA((int) k);
						long elapsed_time_fssa = System.nanoTime() - start_time_fssa;
						System.out.println("elapsed_time_fssa: "+elapsed_time_fssa);
						
						System.out.println("FSSA subset: "+fssa_subset.keySet());
						ArrayList<String> FSSAdominantFeaturesStr = new ArrayList<String>(fssa_subset.keySet());
						ArrayList<Integer> FSSAdominantFeaturesInd = new ArrayList<Integer>();
						int fssa_subset_size = fssa_subset.keySet().size();
						
						for (int i=0; i<FSSAdominantFeaturesStr.size(); i++){ 
							for (int j=0; j<featureNames_Arr_all.length; j++){
								if (FSSAdominantFeaturesStr.get(i).equals(featureNames_Arr_all[j]))
									FSSAdominantFeaturesInd.add(j);
							}
						}
						
						Results res1 = new Results(participant_no, windowSizeSec[w]);
						res1.setAlgNo("1-"+percentages_fssa[perc]*100);
						res1.setWinsizeSec(windowSizeSec[w]);
						res1.setFeaturesNo(fssa_subset_size);
						if (res1.getFeaturesNo() > 0)
							res1.setDominantFeatures(FSSAdominantFeaturesInd);
						res1.setExecTimeFS(elapsed_time_fssa);
						res1.setExecTimeFE(elapsed_time_FE);
						
						results1.add(res1);
					}
					
					
					/**
					 * RELIEF
					 */
					ReliefFAttributeEval relief = new ReliefFAttributeEval();
					Instances ReliefInput = relief.createReliefInput(fssa.getFSSA_input(), featureNames_Arr, finalLabels);
					ArrayList<Integer> ReliefSubset = new ArrayList<Integer>();
					Results res2 = new Results(participant_no, windowSizeSec[w]);
					
					try {
						String[] options = {"-M", "-1", "-D", "1", "-K", "10", "-W", "-A", "2"};
						relief.setOptions(options);
						relief.setSampleSize(-1);
						relief.setSeed(1);
						relief.setNumNeighbours(10);
						long start_time_relief = System.nanoTime();
						relief.buildEvaluator(ReliefInput);
						long elapsed_time_relief = System.nanoTime() - start_time_relief;
						System.out.println("elapsed_time_relief: "+elapsed_time_relief);
						
						ReliefSubset = relief.getIndexes();
						System.out.println("ReliefSubset"+ReliefSubset);
						int Relief_subset_size = ReliefSubset.size();
						res2.setAlgNo("2");
						res2.setWinsizeSec(windowSizeSec[w]);
						res2.setFeaturesNo(Relief_subset_size);
						
						if (!relief.containsNan())
							res2.setDominantFeatures(ReliefSubset);
						else {
							res2.setDominantFeatures(new ArrayList<Integer>());
							res2.setFeaturesNo(0);
						}
						res2.setExecTimeFS(elapsed_time_relief);
						res2.setExecTimeFE(elapsed_time_FE);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
					/**
					 * GCNC
					 */
					System.out.println("GCNC part"+participant_no);
					int nPatterns = finalLabels.size();
					System.out.println("nPatterns: "+nPatterns);	
					
					/**
					 * Initialize GCNC
					 * PREPROCESSING
					 */
					double theta = 0.5;
					Preprocessing prep = new Preprocessing(finalFeatureSet, featureNames_Arr, theta);
					UndirectedGraph graph =  prep.initializeGraph();
	//				p.printGraph();
					prep.writeGraphToFile(graph, "files/features_graph.txt");	
					ArrayList<Integer> graphClusters = prep.louvainClustering();
					System.out.println("graphClusters: "+graphClusters);
					System.out.println(graph.getAllClusters().entrySet());
					
					/**
					 * Start GCNC
					 */
					double delta = 0.0;
					double om = 0;
					GCNC gcnc= new GCNC(finalFeatureSet, featureNames_Arr, nPatterns, graph, delta, om);
					long start_time_GCNC = System.nanoTime();
					gcnc.core();
					long elapsed_time_GCNC = System.nanoTime() - start_time_GCNC;
					System.out.println("elapsed_time_GCNC: "+elapsed_time_GCNC);
					
					ArrayList<Integer> GCNC_dominantFeatInd = gcnc.getReducedSet();
					
					Results res3 = new Results(participant_no, windowSizeSec[w]);
					res3.setAlgNo("3");
					res3.setWinsizeSec(windowSizeSec[w]);
					res3.setFeaturesNo(GCNC_dominantFeatInd.size());
					if (res3.getFeaturesNo() > 0)
						res3.setDominantFeatures(GCNC_dominantFeatInd);
					res3.setExecTimeFS(elapsed_time_GCNC);
					res3.setExecTimeFE(elapsed_time_FE);

					
					results.add(res0);
					results.addAll(results1);
					results.add(res2);
					
					results.add(res3);
					
					
					/**
					 * features file foreach participant
					 */
					Utilities.IO.writeFeatureSetToFile("files/SPL_experiment_per_device/part"+participant_no+"/part"+participant_no+"_dev"+dev+"_winsize"+windowSizeSec[w]+"_features.csv", finalFeatureSet_all,
												finalLabels, featureNames_Arr_all);
				}
			}
			/**
			 * stats file foreach device
			 */
			Utilities.IO.writeResultsToFile("files/SPL_experiment_per_device/stats"+"_dev"+dev+".csv", results);
		}
	}
	
}
