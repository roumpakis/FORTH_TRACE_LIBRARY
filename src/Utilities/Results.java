package Utilities;

import java.util.ArrayList;

public class Results {

	private int participantNo;
	private String AlgNo;
	private int winsizeSec;
	private int featuresNo;
	private ArrayList<Integer> dominantFeatures;
	private long execTimeFS;
	private long execTimeFE;

	public Results(int participantNo, int winsizeSec){
		this.participantNo = participantNo;
		this.winsizeSec = winsizeSec;
		this.dominantFeatures = new ArrayList<Integer>();
	}

	@Override
	public String toString() {
		return "Results [participantNo=" + participantNo + ", AlgNo=" + AlgNo + ", winsizeSec=" + winsizeSec
				+ ", featuresNo=" + featuresNo + ", dominantFeatures=" + dominantFeatures + ", execTimeFS=" + execTimeFS
				+ ", execTimeFE=" + execTimeFE + "]";
	}

	
	public int getParticipantNo() {
		return participantNo;
	}

	public void setParticipantNo(int participantNo) {
		this.participantNo = participantNo;
	}

	public String getAlgNo() {
		return AlgNo;
	}

	public void setAlgNo(String algNo) {
		AlgNo = algNo;
	}

	public int getWinsizeSec() {
		return winsizeSec;
	}

	public void setWinsizeSec(int winsizeSec) {
		this.winsizeSec = winsizeSec;
	}

	public int getFeaturesNo() {
		return featuresNo;
	}

	public void setFeaturesNo(int featuresNo) {
		this.featuresNo = featuresNo;
	}

	public ArrayList<Integer> getDominantFeatures() {
		return dominantFeatures;
	}

	public void setDominantFeatures(ArrayList<Integer> dominantFeatures) {
		this.dominantFeatures = dominantFeatures;
	}

	public long getExecTimeFS() {
		return execTimeFS;
	}

	public void setExecTimeFS(long execTime) {
		this.execTimeFS = execTime;
	}
	
	public long getExecTimeFE() {
		return execTimeFE;
	}

	public void setExecTimeFE(long execTimeFE) {
		this.execTimeFE = execTimeFE;
	}
}
