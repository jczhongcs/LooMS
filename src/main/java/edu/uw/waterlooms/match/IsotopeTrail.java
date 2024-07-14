package edu.uw.waterlooms.match;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;

public class IsotopeTrail implements Comparable<IsotopeTrail>, Serializable {
    public double mz;
    public static double threshold = 0.005; // mz error allowance
    public double maxIntensityRt; // retention time of peak of trail
    public Double[] intensities;
    public Double[] rts;

    public IsotopeTrail() {}

    public IsotopeTrail(double mz, double max, Double[] ints, Double[] rts) {
        this.mz = mz;
        maxIntensityRt = max;
//        this.charge = charge;
        this.intensities = ints;
        this.rts = rts;
    }

    // add rts to this object
    public void ParseRts(String[] rts) {
        this.rts = new Double[rts.length];
        for (int i = 0; i < rts.length; ++i) {
            this.rts[i] = new Double(rts[i]);
        }
    }

    // add intensities to this object
    public void ParseIntensities(String[] intensities) {
        this.intensities = new Double[intensities.length];
        for (int i = 0; i < intensities.length; ++i) {
            this.intensities[i] = new Double(intensities[i]);
        }
    }

    // find peak intensity of trail
    public double FindMaxIntensity() {
        int index = Collections.binarySearch(Arrays.asList(rts), maxIntensityRt);
        return intensities[index];
    }

    @Override
    public int compareTo(IsotopeTrail trail) {
        if (mz > trail.mz + threshold) {
            return 1;
        }
        if (mz < trail.mz - threshold) {
            return -1;
        }
        return 0;
    }

    public String[] peaksToCsv() {
        String[] res = new String[intensities.length];
        for (int i  = 0; i < intensities.length; ++i) {
            res[i] = mz + "," + rts[i] + "," + intensities[i];
        }
        return res;
    }
}
