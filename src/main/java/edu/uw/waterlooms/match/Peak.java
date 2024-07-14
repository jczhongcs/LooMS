package edu.uw.waterlooms.match;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;
import java.util.Comparator;

public class Peak {
    Double rt;
    Double score; // initialized to be peak of score
    Double intensity; // Should be the intensity of the peak
    Double mz;
    String ionType;
    static Comparator<Peak> compareByRt = Comparator.comparing((Peak o) -> o.rt);
    public double getScore(){
        return this.score;
    }
    public void setIonType(String ionType){
        this.ionType = ionType;
    };

    public String getIonType(){
        return this.ionType;
    }

    public Double getRt() {
        return rt;
    }
    public Double getMz() { return mz; }

    public void setIntensity(Double intensity){
        this.intensity = intensity;
    }
    public Double getIntensity() { return this.intensity; };

    /*
    * NOTE: Formerly this was +/- IsotopeTrail.threshold
    * I adjusted this to be +/- 0.5 mz
    * */
    public double compareByMz(Peak o2){
        if (mz > o2.mz + 0.5) {
            return 1;
        }
        if (mz <o2.mz - 0.5) {
            return -1;
        }
        return 0;
    };

    public boolean checkSame(Peak o2){
        if (mz != o2.mz) {
            return false;
        }
        if (rt != o2.rt) {
            return false;
        }
        return true;
    };

    public Peak(double mz, double rt, double score) {
        this.rt = rt;
        this.score = score;
        this.intensity = score;
        this.mz = mz;
    }

    public Peak(double mz, double rt, double score, String ionType) {
        this.rt = rt;
        this.score = score;
        this.intensity = score;
        this.mz = mz;
        this.ionType = ionType;
    }


    // given trail index in a window, construct its peak
    public static Peak getPeakFromIndex(Integer index, IsolationWindow window) {
        XIC xic = window.xics.get(index);
        double maxIntensity = xic.getMaxIntensity();
//        double maxRtRound = Utils.RoundToX(trail.maxIntensityRt, Config.maxGap);
//        double maxRtIntensity = IsolationWindowCollection.max.get(maxRtRound);

        // initialized intensity to be max absolute intensity
        return new Peak(xic.getMZAtMaxIntensity(), xic.getRtAtMaxIntensity(), maxIntensity);
    }

    public String toString() {
        return this.rt + " " + this.mz + " " + this.score;
    }

//    @Override
//    public boolean equals(Object obj) {
//        if (obj == null)
//            return false;
//        Peak other = (Peak)obj;
//        if (Math.abs(mz - other.mz) < IsotopeTrail.threshold && Math.abs(rt - other.rt) < DbMatch.tsThreshold) {
//            return true;
//        }
//        return false;
//    }
}
