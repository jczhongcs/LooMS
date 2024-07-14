package edu.uw.waterlooms.entity;

import edu.uw.waterlooms.match.DbMatch;
import edu.uw.waterlooms.match.Peak;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class IsolationWindowIntegrationTest {

    // TODO: Change this to work on the resource(s) folder
    private static final Double PPM_TOLERANCE = 10d;
    public static final Double mzUpperScalingFactor = 1 + PPM_TOLERANCE * Math.pow(10, -6);
    public static final Double mzLowerScalingFactor = 1 - PPM_TOLERANCE * Math.pow(10, -6);


    // @Test
    public void investigateDuplicateTrails(){
        IsolationWindow isolationWindow = IsolationWindow.parseSerializedIsolationWindowFile(
                "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/serialized_data/",
                1
        );

        DbMatch dbMatch = new DbMatch(null, null, null, null);

        // Analyze WHY this is happening
        double rt = 87.82559;
        double mz = 614.4236;
        double rtTolerance = 01.6;


        assert isolationWindow != null;

        // Sort by RT for subsetting by RT
        isolationWindow.xics.sort(Comparator.comparingDouble(XIC::getRtAtMaxIntensity));
        ArrayList<XIC> prelimTrailSubset = (ArrayList<XIC>) dbMatch.subsetXICBasedOnRT(isolationWindow.xics, rt);


        // Sort by MZ for subsetting by MZ
        prelimTrailSubset.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));

        ArrayList<XIC> finalTrailSubset = (ArrayList<XIC>) dbMatch.subsetXICBasedOnMZ(prelimTrailSubset, mz);
        int b = 0;

        XIC upperXIC = new XIC(mz * mzUpperScalingFactor, rt + rtTolerance);
        XIC lowerXIC = new XIC(mz * mzLowerScalingFactor, rt - rtTolerance);

        int upperIndex = Collections.binarySearch(isolationWindow.xics, upperXIC, Comparator.comparingDouble(XIC::getRtAtMaxIntensity));
        int lowerIndex = Collections.binarySearch(isolationWindow.xics, lowerXIC, Comparator.comparingDouble(XIC::getRtAtMaxIntensity));
        upperIndex = parseInsertionPoint(upperIndex);
        lowerIndex = parseInsertionPoint(lowerIndex);
        ArrayList<XIC> trailsSubset = new ArrayList<> (isolationWindow.xics.subList(lowerIndex, upperIndex));
        int x = 0;


        XIC upperMZXIC = new XIC(mz * mzUpperScalingFactor);
        XIC lowerMZXIC = new XIC(mz * mzLowerScalingFactor);

        //int upperIndex = Collections.binarySearch(isolationWindow.xics, upperMZXIC, Comparator.comparingDouble(XIC::getMZAtMaxIntensity));
        //int lowerIndex = Collections.binarySearch(isolationWindow.xics, lowerMZXIC, Comparator.comparingDouble(XIC::getMZAtMaxIntensity));


        //ArrayList<XIC> trailsSubset = new ArrayList<> (isolationWindow.xics.subList(lowerIndex, upperIndex));

        // subset by MZ first
        Peak upperPeakThreshold = new Peak(mz * mzUpperScalingFactor, 0,0);
        Peak lowerPeakThreshold = new Peak(mz * mzLowerScalingFactor, 0,0);






    }

    private static int parseInsertionPoint(int index){
        if (index < 0){
            return (-1 * (index) - 1);
        } else
        {
            return index;
        }
    }
}
