package Utilities;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

public class IO {

	public static void writeFeatureSetToFile(String filename, HashMap<String, ArrayList<Double>> features,
			ArrayList<Double> labels, String[] featureNames) {
		
		try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename, true)))) {
			
			int sz = features.get(featureNames[0]).size();
			
			for (int i=0; i<featureNames.length; i++){
				out.print(featureNames[i]+",");
			}
			out.println("label");
			for (int i=0; i<featureNames.length; i++){
				out.print(i+",");
			}
			out.println(000);
			
			for (int i=0; i<sz; i++){
				for (int j=0; j<features.size(); j++){
					out.print(features.get(featureNames[j]).get(i)+",");
				}
				out.println(labels.get(i));
			}

		} catch (IOException e) {
			// exception handling left as an exercise for the reader
		}
	}
	public static void writeResultsToFile(String filename, ArrayList<Results> results){
		try {
			
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(filename, true)));
			out.println("participant,alg,winsize,featuresNo,dominantFeatures,execTimeFSA,exectimeFE");
			for (int i=0; i<results.size(); i++){
				Results res = results.get(i);
				out.print(res.getParticipantNo()+",");
				out.print(res.getAlgNo()+",");
				out.print(res.getWinsizeSec()+",");
				out.print(res.getFeaturesNo()+",");
				ArrayList<Integer> domFeats = res.getDominantFeatures();
				for (int k=0; k<domFeats.size(); k++){
					out.print(domFeats.get(k)+",");
				}
				out.print(res.getExecTimeFS()+",");
				out.print(res.getExecTimeFE());
				out.println();
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
