package edu.uw.waterlooms.match;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import static java.lang.Math.log10;
import static java.lang.Math.max;

public class DbMatch implements Comparable<DbMatch>{

    private static final double PPM_TOLERANCE = 10d;
    public static final double mzUpperScalingFactor = 1 + PPM_TOLERANCE * Math.pow(10, -6);
    public static final double mzLowerScalingFactor = 1 - PPM_TOLERANCE * Math.pow(10, -6);

    public static final double RT_TOLERANCE_MINUTES = 0.16;

    public MatchedInterval matchedRts; // stores the interval where peptide most likely occurs
    public String id; // id of protein
    public String composition; // peptide composition
    private final Peptide peptide; // peptide object that is matched

    public DbMatch(String id, String composition, MatchedInterval matchedRts, Peptide p) {
        this.composition = composition;
        this.matchedRts = matchedRts;
        this.id = id;
        this.peptide = p;
    }

    public Peptide getPeptide(){
        return this.peptide;
    }

    @Override
    public int compareTo(DbMatch other) {
        Double score = this.matchedRts.score;
        return score.compareTo(other.matchedRts.score);
    }

    public static double calculateScoreContribution(double x) {
        // Implement score f'n from class
        double scoreContribution = 0d;
        if (x > 0.01){
            scoreContribution = log10(100 * x);
        }
        return max(0d, scoreContribution);
    }

    /**
     * Given an array of of mz's...
     * Given an Isolation Window with TRAILs ...
     * Subset the TRAILs by each mz value in the array.
     * TODO :THIS IS INCREDIBLY UNOPTIMIZED
     *
     * @param mzArray double[] of mz values to select TRAILs with
     * @param trails ArrayList of Trails to subset from
     * @return ArrayList<Peak>
     */
    public static ArrayList<XIC> subsetXICbyMZArray(
            double[] mzArray,
            ArrayList<XIC> trails
    ){
        ArrayList<XIC> extractedXIC = new ArrayList<>();

        // For each mz value in the array, subset the XIC based on the mz value alone, return a list of PEAK representations of each matched XIC
        for (double mz : mzArray) {
            List<XIC> mzExtractedTrails = subsetXICBasedOnMZ(trails, mz);
            // TODO: To cut out on the MatchFeatures.matchSeries here, you can just do the following:
            // Take the max of the mzExtracted Trails
            // Avoid the Collections.addAll call altogether
            // For the given label & mz, add it to the MatchFeatures object


            //  runtime 4m 27s 486ms
            Collections.addAll(extractedXIC, mzExtractedTrails.toArray(new XIC[0]));

            //  runtime 4m 24s 856ms
            // extractedXIC.addAll(mzExtractedTrails);


            // runtime (with sort) 4m 54s 79ms
        }
        return extractedXIC;
    }


    public static DbMatch SPScoredMatchedIntervals(
            ArrayList<XIC> extractedTrails,
            Peptide peptide,
            IsolationWindow isolationWindow,
            ArrayList<XIC> precursors,
            boolean outputDebug,
            String debugOutfile
    ){
        ArrayList<MatchedInterval> allMatchedIntervals = new ArrayList<>();

        // Compute set of Precursors at charge 2
        double mz = Utils.MassToMz(peptide.mass, 2);
        List<XIC> filteredPrecursors = subsetXICBasedOnMZ(precursors, mz);

        // For a given Peptide, and its extracted peaks ...
        for (XIC trail: extractedTrails){
            List<XIC> peaksWithinRTWindow = subsetXICBasedOnRT(extractedTrails, trail.getRtAtMaxIntensity());

            // Compute the "SP" score for the MatchedInterval
            double spscore = 0;
            spscore += computeSPForIonSeries((ArrayList<XIC>) peaksWithinRTWindow, peptide.y_ions, isolationWindow, 'y');
            spscore += computeSPForIonSeries((ArrayList<XIC>) peaksWithinRTWindow, peptide.b_ions, isolationWindow, 'b');

            MatchedInterval matchedInterval = null;

            // This is is precursor matching in every window
            ArrayList<XIC> matchedPrecursors = binarySearchPrecursorList(filteredPrecursors, mz, trail.getRtAtMaxIntensity());
            // if matchedPrecursors is not empty
            if (matchedPrecursors.size() > 0){
                // maxIntensity is already normalized to intensity_window_evg in the reading f'n
                // Take the log10 of it
                spscore += calculateScoreContribution(matchedPrecursors.get(0).getMaxIntensity());
                matchedInterval = new MatchedInterval(spscore, peaksWithinRTWindow, matchedPrecursors.get(0));
            } else {
                matchedInterval = new MatchedInterval(spscore, peaksWithinRTWindow);
            }
            allMatchedIntervals.add(matchedInterval);
        }

        if(outputDebug){
            // Call the write helper f'n here
            try {
            writeMatchWindows(allMatchedIntervals, debugOutfile);
                } catch (IOException $e){
                System.err.println($e.getMessage());
            }
        }

        // Stream to get the max interval in set
        MatchedInterval maxScoreInterval = Collections.max(allMatchedIntervals, Comparator.comparing(x -> x.score));
        return new DbMatch(peptide.id, peptide.composition, maxScoreInterval, peptide);
    }

    public static void writeMatchWindows(
            ArrayList<MatchedInterval> matchedIntervals,
            String outputFile
    ) throws IOException {
        // Setup the Output file
        FileWriter fileWriter = new FileWriter(outputFile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
        // Write the Header
        // rtAtFirstPeak, score,
        bufferedWriter.write(
                "rtAtFirstPeak" + "\t" +
                        "logIntensityScore" + "\t" +
                        "score" + "\t" +
                        "matchedPrecursor" + "\n"
        );

        for (MatchedInterval match : matchedIntervals){
            if(match.matchedFeatures == null){
                continue;
            }
            bufferedWriter.write(match.matchedFeatures.allRetentionTimes.get(0) + "\t");
            bufferedWriter.write(match.matchedFeatures.maxPeakLogIntensityScore + "\t");
            bufferedWriter.write(match.matchedFeatures.score + "\t");
            if(match.precursor != null){
                bufferedWriter.write("TRUE");
            } else {
                bufferedWriter.write("FALSE");
            }
            bufferedWriter.write("\n");
        }
    }


    private static ArrayList<XIC> binarySearchPrecursorList(List<XIC> precursors, double mz, double rt){
        final double PPM_TOLERANCE = 30;
        XIC upperMZXIC = new XIC(mz * (1 + PPM_TOLERANCE * Math.pow(10, -6)));
        XIC lowerMZXIC = new XIC(mz * (1 - PPM_TOLERANCE * Math.pow(10, -6)));

        // Heuristically cut the search space down just based on ppm
        int upperIndex = Collections.binarySearch(precursors, upperMZXIC, XIC::compareMzXIC);
        int lowerIndex = Collections.binarySearch(precursors, lowerMZXIC, XIC::compareMzXIC);
        upperIndex = parseInsertionPoint(upperIndex);
        lowerIndex = parseInsertionPoint(lowerIndex);

        double upper = rt + Config.tsThreshold;
        double lower = rt - Config.tsThreshold;

        List<XIC> possiblePrecursors = precursors.subList(lowerIndex, upperIndex);

        ArrayList<XIC> matchedPrecursors = possiblePrecursors.stream().filter(
                xic -> (xic.getRtAtMaxIntensity() <= upper) && (lower <= xic.getRtAtMaxIntensity())
        ).collect(Collectors.toCollection(ArrayList::new));

        return matchedPrecursors;
    }

    /**
     * Given an ArrayList of XIC, subset the ArrayList by RT and some tolerance defined in the constant.
     * @param trails ArrayList of XIC objects to be subsetted
     * @param rt double RetentionTime to subset against
     * @return List of XIC objects that are extracted from the superset trails
     */
    public static List<XIC> subsetXICBasedOnRT(ArrayList<XIC> trails, double rt){
        XIC upperRTXIC = new XIC(0d, rt + RT_TOLERANCE_MINUTES);
        XIC lowerRTXIC = new XIC(0d, rt - RT_TOLERANCE_MINUTES);

        int upperIndex = Collections.binarySearch(trails, upperRTXIC, XIC::compareRtXIC);
        int lowerIndex = Collections.binarySearch(trails, lowerRTXIC, XIC::compareRtXIC);
        upperIndex = parseInsertionPoint(upperIndex);
        lowerIndex = parseInsertionPoint(lowerIndex);

        return new ArrayList<> (trails.subList(lowerIndex, upperIndex));
    }

    /**
     * Given an ArrayList of XIC, subset the ArrayList by m/z and some ppm value defined in the constant.
     * @param trails ArrayList of XIC objects to be subsetted
     * @param mz double m/z value to subset against
     * @return List of XIC objects that are extracted from the superset trails
     */
    public static List<XIC> subsetXICBasedOnMZ(ArrayList<XIC> trails, Double mz){
        XIC upperMZXIC = new XIC(mz * mzUpperScalingFactor);
        XIC lowerMZXIC = new XIC(mz * mzLowerScalingFactor);

        int upperIndex = Collections.binarySearch(trails, upperMZXIC, XIC::compareMzXIC);
        int lowerIndex = Collections.binarySearch(trails, lowerMZXIC, XIC::compareMzXIC);
        upperIndex = parseInsertionPoint(upperIndex);
        lowerIndex = parseInsertionPoint(lowerIndex);

        return new ArrayList<> (trails.subList(lowerIndex, upperIndex));
    }

    /**
     * Process the (possibly negative) index returned by Collections.binarySearch.
     * @param index int index of the binary searched item.
     * @return the index of the searched item OR the index to insert the new item.
     */
    private static int parseInsertionPoint(int index){
        if (index < 0){
            return (-1 * (index) - 1);
        } else
        {
            return index;
        }
    }

    /**
     * Match a peptide against the superset of Trails.
     * @param peptide Peptide to be matched.
     * @param isolationWindow IsolationWindow that the peptide falls into by precursor at z=1.
     * @param precursors ArrayList<XIC> precursors in this window
     * @return
     */
    public static ArrayList<MatchedInterval> peptideMatch(
            Peptide peptide,
            IsolationWindow isolationWindow,
            IsolationWindow precursorIsolationWindow,
            ArrayList<XIC> precursors
    ){
        String p = peptide.composition;
        if (p.equals("VTNPNAR")) {
            int jjj=0;
        }
        ArrayList<MatchedInterval> allIntervals = new ArrayList<>();
        allIntervals.ensureCapacity(isolationWindow.rtToTrailSet.size() + 1); // Pre-load the size for optimization

        isolationWindow.rtToTrailSet.forEach((rt, trails) -> {
            // Raw-Score the peptide matchFeatures
            MatchFeatures matchFeatures = new MatchFeatures(peptide, 1);
            matchFeatures.matchIons(trails, isolationWindow.trailToMatchedPeptidesSet);
            // matchPeptide here from the trails given the RT
            matchFeatures.matchPrecursor(isolationWindow.rtToPrecursorSet.get(rt));

            matchFeatures.scoreMatchInterval(isolationWindow, precursorIsolationWindow, trails);

            if (matchFeatures.precursor != null) {
                MatchedInterval interval = new MatchedInterval(matchFeatures);
                interval.setScore(matchFeatures.score);
                interval.setPrecursor(matchFeatures.precursor);

                // Append the basic interval to the set of intervals
                allIntervals.add(interval);
            }

            // 1m 1s for 500 peptides

            /*
            // Given the peptide's b & y ions, match them to the trails within the RT window
            ArrayList<XIC> allTrails = subsetXICbyMZArray(peptide.y_ions, trails);

            // 2m 11s 769ms runtime
            Collections.addAll(allTrails, subsetXICbyMZArray(peptide.b_ions, trails).toArray(new XIC[0]));

            // 2m 9s 480ms runtime
            //allTrails.addAll(subsetXICbyMZArray(peptide.b_ions, trails));


            if (allTrails.size() != 0){
                // Raw-Score the peptide matchFeatures
                MatchFeatures matchFeatures = new MatchFeatures(peptide, 1);
                matchFeatures.matchIons(allTrails);
                matchFeatures.scoreMatchInterval(isolationWindow, allTrails);
                MatchedInterval interval = new MatchedInterval(matchFeatures);
                interval.setScore(matchFeatures.score);

                // Append the basic interval to the set of intervals
                allIntervals.add(interval);
            }

             */

        });

        return allIntervals;
    }

    private static double computeSPForIonSeries(
            ArrayList<XIC> ions,
            double[] ion_mz,
            IsolationWindow isolationWindow,
            char type
    ) {
        boolean consecutive = false;
        double score = 0;
        double reward = 0;
        int m = 0;

        ArrayList<Double> allRTs = new ArrayList<Double>();

        for (double mz : ion_mz){
            if (mz == 0){
                // Skip the 0 value
                continue;
            }

            List<XIC> match = subsetXICBasedOnMZ(ions, mz);

            if (!match.isEmpty()){
                // If the resultant matches for the particular fragment ion is non null
                // Append the intensity to the finalized b_score
                // formerly .stream().max
                //XIC maxIntensityPeak = match.stream().max(Comparator.comparing(XIC::getMaxIntensity)).get();

                // Updated to Collections.max() to avoid the .stream().max call
                XIC maxIntensityPeak = Collections.max(match, Comparator.comparingDouble(XIC::getMaxIntensity));


                double rt = maxIntensityPeak.getRtAtMaxIntensity();
                allRTs.add(rt);
                double intensity = maxIntensityPeak.getMaxIntensity();

                // Compute the RANK of a given peak
                ArrayList<Double> intensityRanks = isolationWindow.rtToIntensitySet.get(rt);
                double rank = Collections.binarySearch(intensityRanks, intensity);
                double rankValue = rank / intensityRanks.size();
                score += rankValue;

                // Search for a neutral loss here
                double waterLossMz = computeWaterNeutralLoss(mz);
                double ammoniaLossMz = computeAmmoniaNeutralLoss(mz);
                double isotopicPlusOneMz = computeIsotopicMz(mz);

                // Search for neutral losses
                // FEATURE b/y - 18
                List<XIC> waterMatch = subsetXICBasedOnMZ(ions, waterLossMz);
                if(!waterMatch.isEmpty()) {
                    // Naively get 1st peak
                    XIC waterLossPeak = waterMatch.get(0);
                    double normalizedWaterlossIntensity = normalizeAgainstMaxIntensity(waterLossPeak, isolationWindow);
                    score += calculateScoreContribution(normalizedWaterlossIntensity);
                }
                // FEATURE b/y - 17
                List<XIC> ammoniaMatch = subsetXICBasedOnMZ(ions, ammoniaLossMz);
                if(!ammoniaMatch.isEmpty()) {
                    // Naively get 1st peak
                    XIC ammoniaLossPeaks = ammoniaMatch.get(0);
                    double normalizedAmmonialossIntensity = normalizeAgainstMaxIntensity(ammoniaLossPeaks, isolationWindow);
                    score += calculateScoreContribution(normalizedAmmonialossIntensity);
                }


                // FEATURE: Normalized Intensity
                double normalizedIntensity = normalizeAgainstMaxIntensity(maxIntensityPeak, isolationWindow);
                score += calculateScoreContribution(normalizedIntensity);
                consecutive = true;
                m += 1;
            } else {
                consecutive = false;
            }
            if (consecutive){
                reward += 0.075;
            }
        }

        double intermediateScore = (score * m * (1 + reward)) / (ion_mz.length - 1);

        // If there are matches, help to add some score to the interval
        if (allRTs.size() > 0){
            double maxRT = Collections.max(allRTs);
            double minRT = Collections.min(allRTs);
            intermediateScore += (
                    ((2 * RT_TOLERANCE_MINUTES) - (maxRT - minRT)) / (2 * RT_TOLERANCE_MINUTES)
            );
        }

        // Formerly 2* ion_mz.length, but thats fine if its not, don't want to unfairly penalize LONG peptide(s)
        return intermediateScore;
    }

    private static double normalizeAgainstMaxIntensity(XIC trail, IsolationWindow isolationWindow){
        double rt = trail.getRtAtMaxIntensity();
        double intensity = trail.getMaxIntensity();
        double intensityNormalization = isolationWindow.max.get(rt);

        return (intensity / intensityNormalization);
    }

    private static double computeWaterNeutralLoss(double mz){
        // Implicitly only works for m/z at z = 1
        double waterMass = 18.0153;
        return (mz - waterMass);
    }
    private static double computeAmmoniaNeutralLoss(double mz){
        // Implicitly only works for m/z at z = 1
        double ammoniaMass = 17.02654;
        return (mz - ammoniaMass);
    }

    private static double computeIsotopicMz(double mz){
        // Implicitly only works for m/z at z = 1
        double daltonMass = 1.00727647;
        return (mz + daltonMass);
    }

}
