package edu.uw.waterlooms.peptideMatch;


import java.util.Comparator;

public class Peak {
    Double rt;
    Double score; // initialized to be peak of score
    Double mz;
    static Comparator<Peak> compareByRt = Comparator.comparing((Peak o) -> o.rt);
//    static Comparator<Peak> compareByRtMz = (Peak o1, Peak o2) -> o1.compareByRtMzFn(o2);
//
//    public int compareByRtMzFn(Peak o2){
//        int rtRes = rt.compareTo(o2.rt);
//        if (rtRes != 0) {
//            return rtRes;
//        }
//        return mz.compareTo(o2.mz);
//    };

    public double compareByMz(Peak o2){
        if (mz > o2.mz + MSOneTrail.threshold) {
            return 1;
        }
        if (mz <o2.mz - MSOneTrail.threshold) {
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
        this.mz = mz;
    }

    // given trail index in a window, construct its peak
/*
    public static Peak getPeakFromIndex(Integer index, MSTwoTrailSet window) {
        MSTwoTrail trail = window.arrMSTwoTrail.get(index);
        double maxIntensity = trail.FindMaxIntensity();
//        double maxRtRound = Utils.RoundToX(trail.maxIntensityRt, Config.maxGap);
//        double maxRtIntensity = IsolationWindowCollection.max.get(maxRtRound);

        // initialized intensity to be max absolute intensity

        return new Peak(trail.mz, trail.maxIntensityRt, maxIntensity);
    }

    public String toString() {
        return this.rt + " " + this.mz + " " + this.score;
    }
*/
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
