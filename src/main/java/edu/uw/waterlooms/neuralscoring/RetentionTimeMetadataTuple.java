package edu.uw.waterlooms.neuralscoring;

import java.io.Serializable;
import java.util.ArrayList;

public class RetentionTimeMetadataTuple  implements Serializable {
    /*
    * Maximum Matched Peak Intensity
    * Maximum UnMatched Peak Intensity
    * Total Ion Current
    * Number Matched Ions
    * */
    public double maximumMatchedPeakIntensity;
    public double maximumUnmatchedPeakIntensity;
    public double totalIonCurrent;
    public int numberMatchedIons;

    public double[] maximumMatchedPeakIntensities;
    public double[] maximumUnmatchedPeakIntensities;

    /**
     * Lightweight constructor
     */
    public RetentionTimeMetadataTuple(){
        this.maximumMatchedPeakIntensity = 0d;
        this.maximumUnmatchedPeakIntensity = 0d;
        this.totalIonCurrent = 0d;
        this.numberMatchedIons = 0;

        this.maximumMatchedPeakIntensities = new double[2000];
        this.maximumUnmatchedPeakIntensities = new double[2000];
    }

    public void setMaximumMatchedPeakIntensity(double maximumMatchedPeakIntensity) {
        this.maximumMatchedPeakIntensity = maximumMatchedPeakIntensity;
    }

    public void setMaximumUnmatchedPeakIntensity(double maximumUnmatchedPeakIntensity){
        this.maximumUnmatchedPeakIntensity = maximumUnmatchedPeakIntensity;
    }

    public void setTotalIonCurrent(double totalIonCurrent){
        this.totalIonCurrent = totalIonCurrent;
    }

    public void setNumberMatchedIons(int numberMatchedIons){
        this.numberMatchedIons = numberMatchedIons;
    }

    public double getMaximumMatchedPeakIntensity(){
        return this.maximumMatchedPeakIntensity;
    }

    public double getMaximumUnmatchedPeakIntensity(){
        return this.maximumUnmatchedPeakIntensity;
    }

    public double getTotalIonCurrent(){
        return this.totalIonCurrent;
    }

    public int getNumberMatchedIons(){
        return this.numberMatchedIons;
    }

    // TODO: Indexed by 0 so need to figure this out
   public void updateMaxMatchedGivenMZ(double mz, double intensity){
       // Get the floor of the mz value
       // Have a statement that adjusts it to the 2000 bin if m/z > 2000
       // Check the intensity value in the bin
       // Update if it's larger
       int bin = (int) Math.floor(mz);
       bin = (bin >= 2000) ? 1999 : bin; // Sanity check to fix the bin value here

       // Update the intensity value at the bin if it exists
       if (this.maximumMatchedPeakIntensities[bin] < intensity){
           this.maximumMatchedPeakIntensities[bin] = intensity;
       }
    }

    public void updateMaxUnmatchedGivenMZ(double mz, double intensity){
        // Get the floor of the mz value
        // Have a statement that adjusts it to the 2000 bin if m/z > 2000
        // Check the intensity value in the bin
        // Update if it's larger
        int bin = (int) Math.floor(mz);
        bin = (bin >= 2000) ? 1999 : bin; // Sanity check to fix the bin value here
        // Update the intensity value at the bin if it exists
        if (this.maximumUnmatchedPeakIntensities[bin] < intensity){
            this.maximumUnmatchedPeakIntensities[bin] = intensity;
        }
    }

    public double[] getMaximumMatchedPeakIntensities() {
        return maximumMatchedPeakIntensities;
    }

    public double[] getMaximumUnmatchedPeakIntensities() {
        return maximumUnmatchedPeakIntensities;
    }
}