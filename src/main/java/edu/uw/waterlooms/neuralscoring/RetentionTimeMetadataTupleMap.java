package edu.uw.waterlooms.neuralscoring;

import edu.uw.waterlooms.entity.XIC;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class RetentionTimeMetadataTupleMap {
    public HashMap<Double, RetentionTimeMetadataTuple> retentionTimeMetadataTupleHashMap;
    public RetentionTimeMetadataTupleMap(){
        this.retentionTimeMetadataTupleHashMap = new HashMap<>();
    }
    public void addRetentionTimeToMetadataElement(double retentionTime, RetentionTimeMetadataTuple tuple){
        this.retentionTimeMetadataTupleHashMap.put(retentionTime, tuple);
    }
    public boolean retentionTimeExists(double retentionTime){
        return this.retentionTimeMetadataTupleHashMap.containsKey(retentionTime);
    }

    public void updateSetGivenTrail(XIC trail, int count){
        ArrayList<Double> retentionTimes = trail.getRetentionTimes();
        ArrayList<Double> intensities = trail.getIntensities();
        ArrayList<Double> mzs = trail.getMassChargeRatios();

        for(int trailIdx = 0; trailIdx < retentionTimes.size(); trailIdx ++){
            double rt = retentionTimes.get(trailIdx);
            double intensity = intensities.get(trailIdx);
            double mz = mzs.get(trailIdx);


            // If RT does NOT exist, create it
            if (!retentionTimeExists(rt)){
                // Create the RetentionTimeMetadataTuple
                RetentionTimeMetadataTuple newRtTuple = new RetentionTimeMetadataTuple();
                if (count > 0){
                    newRtTuple.numberMatchedIons = count;
                    newRtTuple.maximumMatchedPeakIntensity = intensity;
                    newRtTuple.totalIonCurrent = intensity;
                } else {
                    newRtTuple.numberMatchedIons = 0;
                    newRtTuple.maximumUnmatchedPeakIntensity = intensity;
                    newRtTuple.totalIonCurrent = intensity;
                }
                addRetentionTimeToMetadataElement(rt, newRtTuple);
            } else {
                // Update the Existing Metadatatuple
                RetentionTimeMetadataTuple retrievedRtTuple = this.retentionTimeMetadataTupleHashMap.get(rt);
                if (count > 0){
                    // UpdateMaximumMatchedPeakIntensity
                    retrievedRtTuple.setNumberMatchedIons(retrievedRtTuple.getNumberMatchedIons() + count); // +1 to num matches
                    if (intensity > retrievedRtTuple.getMaximumMatchedPeakIntensity()){
                        retrievedRtTuple.setMaximumMatchedPeakIntensity(intensity);
                    }
                    retrievedRtTuple.setTotalIonCurrent(retrievedRtTuple.getTotalIonCurrent() + intensity);
                    retrievedRtTuple.updateMaxMatchedGivenMZ(mz, intensity);
                } else {
                    // Update MaximumUnmatchedPeakIntensity
                    if (intensity > retrievedRtTuple.getMaximumUnmatchedPeakIntensity()){
                        retrievedRtTuple.setMaximumUnmatchedPeakIntensity(intensity);
                    }
                    retrievedRtTuple.setTotalIonCurrent(retrievedRtTuple.getTotalIonCurrent() + intensity);
                    retrievedRtTuple.updateMaxUnmatchedGivenMZ(mz, intensity);
                }
            }
        }
    }

    public void normalizeRtTupleMapElements(HashMap<Double, Double> rtToMaximumIntensityMap){
        for (HashMap.Entry<Double, Double> entry  : rtToMaximumIntensityMap.entrySet()) {
            Double rt = entry.getKey();
            Double maximumIntensityGivenRT = entry.getValue();

            RetentionTimeMetadataTuple tuple = this.retentionTimeMetadataTupleHashMap.get(rt);

            // sanity check
            if (tuple == null){
                continue;
            }

            // For each value in max matched & unmatched sets, normalize to max
            for (int i = 0; i < tuple.getMaximumMatchedPeakIntensities().length; i ++){
                tuple.getMaximumMatchedPeakIntensities()[i] /= maximumIntensityGivenRT;
                tuple.getMaximumUnmatchedPeakIntensities()[i] /= maximumIntensityGivenRT;
            }
        }
    }

    public void writeMapToTSV(int windowIdx, String workingDirectory) throws IOException {
        FileWriter fileWriter = new FileWriter(workingDirectory + windowIdx + "_retentionTimeMzIntensityBinMap.tsv");
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);


        bufferedWriter.write("rt" + "\t");
        for (int i = 0; i < 2000; i++){
            bufferedWriter.write(i + "_da_max_matched" + "\t");
        }
        for (int i = 0; i < 2000; i++){
            bufferedWriter.write(i + "_da_max_unmatched" + "\t");
        }
        bufferedWriter.write("\n");


        // For each RT and MetadataTuple in the map
        this.retentionTimeMetadataTupleHashMap.forEach((rt, metadata) ->{
            try {
                bufferedWriter.write(Double.toString(rt) + "\t");
                double[] matched = metadata.getMaximumMatchedPeakIntensities();
                double[] unmatched = metadata.getMaximumUnmatchedPeakIntensities();
                for (int i = 0; i< matched.length; i ++){
                    bufferedWriter.write(Double.toString(matched[i]) + "\t");
                }
                for (int i = 0; i< unmatched.length; i ++){
                    bufferedWriter.write(Double.toString(unmatched[i]) + "\t");
                }
                bufferedWriter.write("\n");
                //bufferedWriter.write(Double.toString(metadata.getMaximumMatchedPeakIntensity()) + "\t");
                //bufferedWriter.write(Double.toString(metadata.getMaximumUnmatchedPeakIntensity()) + "\t");
                //bufferedWriter.write(Double.toString(metadata.getTotalIonCurrent()) + "\t");
                //bufferedWriter.write(Integer.toString(metadata.numberMatchedIons) + "\n");
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
        bufferedWriter.close();
        fileWriter.close();
    }


}
