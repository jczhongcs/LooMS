package edu.uw.waterlooms.utility;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.Pair;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.match.*;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;

public class GeneratePeptideTrailMatchAcrossAllRTsIntegrationTest {
    //final String PEPTIDE_OF_INTEREST = "LDVGNAEVK"; First peptide at <1% FDR
    // final String PEPTIDE_OF_INTEREST = "FDSVGGLSNHIAALK"; // Arbitrary peptide
    final String PEPTIDE_OF_INTEREST ="VGGNIEVLGFNAR"; // Peptide with bad trail match at MSGF+ reported RT (r01,q01)
    final String SERIALIZED_DATA_DIR = "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/output/r01_trails/";
    final String PRECURSOR_FILE = "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/r01_dia_data.precursors";
    final String ISOLATION_WINDOW_RANGE_FILE = SERIALIZED_DATA_DIR + "isolationWindowRanges.out";
    final String DEBUG_DIR = "/home/jia/Documents/Code/phd_research/data/trail_debug/";

//    @Test
    public void generateXICsForOnePeptideOfInterest() throws IOException{
        Peptide peptide = new Peptide("", this.PEPTIDE_OF_INTEREST);
        ArrayList<XIC> precursors = RunFilter.parsePrecursorDataset(
                PRECURSOR_FILE,
                2
        );

        int windowIndex = 1;
        RunFilter runFilter = new RunFilter();
        ArrayList<Pair> windowRanges= IsolationWindow.parseWindowRanges(ISOLATION_WINDOW_RANGE_FILE);

        for (Pair window : windowRanges){
            //
            double tmpMz = Utils.MassToMz(peptide.getMass(), 2);
            double lBound = new Double(window.getL().toString());
            double rBound = new Double(window.getR().toString());
            if (lBound <= tmpMz && tmpMz <= rBound){
                break;
            }
            windowIndex += 1;
        }
        IsolationWindow deserializedWindowAndTrails = IsolationWindow.parseSerializedIsolationWindowFile(this.SERIALIZED_DATA_DIR, windowIndex);
        assert deserializedWindowAndTrails != null;
        deserializedWindowAndTrails.xics.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));

        ArrayList<MatchedInterval> matchedIntervalsForPeptide =
                DbMatch.peptideMatch(peptide, deserializedWindowAndTrails, deserializedWindowAndTrails, precursors); // TODO THIS IS A BUG the 2nd deserializedWindowAndTrails should not be this

        // Sort the MatchedIntervalsForPeptide since they might be keyed differently
        matchedIntervalsForPeptide.sort(Comparator.comparingDouble(
                matchedInterval -> {
                    if(matchedInterval.getMatchedFeatures().allRetentionTimes.size() != 0){
                        return matchedInterval.getMatchedFeatures().allRetentionTimes.get(0);
                    } else {
                        return 0;
                    }
                })
        );

        // Write all the matchedIntervals
        // Write a pseudo _result.tsv file here
        runFilter.writeAllMatchedIntervalsAndDummyResultTSVForAGivenPeptide(
                matchedIntervalsForPeptide,
                DEBUG_DIR,
                peptide.composition
        );
        runFilter.writeMatchIntervalMatchedFeaturesTSV(matchedIntervalsForPeptide,DEBUG_DIR);



    }

}
