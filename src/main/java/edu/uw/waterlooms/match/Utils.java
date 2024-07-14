package edu.uw.waterlooms.match;

import static java.lang.Math.round;

public class Utils {
    public static int DeterminePepCharge(double mass) {
        if (mass <= 800) {
            return 2;
        }
        return 1;
    }

    // x is how much you want to round to (1/3 etc)
    public static double RoundToX(double o, double x) {
        if (x == 0) {
            return o;
        }
        double res = round(o / x) * x;
        return res;
    }


    public static double MassToMz(double mass, int z) {
        return mass/z + AminoAcid.ProtonMass;
    }

    public double MzToMass(double mz, int z) {
        return (mz - AminoAcid.ProtonMass) * z;
    }

}
