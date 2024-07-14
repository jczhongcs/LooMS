package edu.uw.waterlooms.peptideMatch;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.log10;

public abstract class DbBaseMatch {
    public String id; // id of protein
    public String composition; // peptide composition
    public List<Integer> bpos = new ArrayList<>();
    public List<Integer> ypos = new ArrayList<>();

    public List<Integer> bipos = new ArrayList<>();
    public List<Integer> yipos = new ArrayList<>();

    public double score;


    public static double EmpiricalScore(double intensity) {
        // so score is >=0
        if (intensity < 1) {
            return 0;
        }
        return log10(intensity);

//        return 1;
    }


//    public abstract DbBaseMatch  PepMatch(Peptide pep, List<MSTwoTrail> lMSTwotrails,double maxPeakArea);
    public abstract MatchPeptdide PepMatch(Peptide pep, List<MSTwoTrail> lMSTwotrails,
                                           double maxPeakArea, Enums.ScoreMode scoreMode, BufferedWriter writer, int icount,
                                           double ms1rt, double ms1Mass, long ms1id,double ms1qs) ;


    protected static int binarySearch0(List<MSTwoTrail>  a, int fromIndex, int toIndex,
                                     double key) {
        int low = fromIndex;
        int high = toIndex - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            double midVal = a.get(mid).getMzApex();

            if (midVal < key)
                low = mid + 1;  // Neither val is NaN, thisVal is smaller
            else if (midVal > key)
                high = mid - 1; // Neither val is NaN, thisVal is larger
            else {
                long midBits = Double.doubleToLongBits(midVal);
                long keyBits = Double.doubleToLongBits(key);
                if (midBits == keyBits)     // Values are equal
                    return mid;             // Key found
                else if (midBits < keyBits) // (-0.0, 0.0) or (!NaN, NaN)
                    low = mid + 1;
                else                        // (0.0, -0.0) or (NaN, !NaN)
                    high = mid - 1;
            }
        }

        return -(low + 1);  // key not found.
    }

}
