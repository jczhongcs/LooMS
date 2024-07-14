package edu.uw.waterlooms.match;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.service.PeptideScoreService;

import java.util.*;
import java.util.stream.Collectors;

public class MatchFeatures {
  private static final double PPM_TOLERANCE = 10d;
  public static final double mzUpperScalingFactor = 1 + PPM_TOLERANCE * Math.pow(10, -6);
  public static final double mzLowerScalingFactor = 1 - PPM_TOLERANCE * Math.pow(10, -6);
  private final PeptideScoreService scoreService;

  public double score;
  public double neutralLossScore;
  public double precursorScore;
  public double rtDistributionScore;
  public double spScore;
  public double rankScore;
  public double maxPeakLogIntensityScore;
  public double trailLogIntensityScore;
  public double peaksPerTrailScore;
  public double rtIntervalIntensityScore;
  public int numIonsMatched;
  public double rtStandardDeviation;

  public Peptide peptide;
  public XIC precursor;
  protected int fragmentIonCharge;
  public ArrayList<Double> allRetentionTimes; // For computing average could possibly remove

  // All matches are done from the N->C direction
  public ArrayList<XIC> bIonsMatched;
  public ArrayList<XIC> yIonsMatched;
  public ArrayList<XIC> bIonsPlusOneMatched; // +C13 Isotope Match
  public ArrayList<XIC> yIonsPlusOneMatched; // +C13 Isotope Match
  public ArrayList<XIC> bIonH2OLossMatched;
  public ArrayList<XIC> bIonNH3LossMatched;
  public ArrayList<XIC> yIonH2OLossMatched;
  public ArrayList<XIC> yIonNH3LossMatched;

  public MatchFeatures(Peptide peptide, int fragmentIonCharge) {
    this.peptide = peptide;
    this.fragmentIonCharge = fragmentIonCharge;
    this.scoreService = new PeptideScoreService();
    this.allRetentionTimes = new ArrayList<>();
    this.allRetentionTimes.ensureCapacity(50); // Todo pre-set at 50 to reduce the ensureCapacityInternal call
    this.numIonsMatched = 0;
  }


  public void scoreMatchInterval(IsolationWindow isolationWindow, IsolationWindow precursorIsolationWindow,List<XIC> trailsWithinRTWindow) {
    // For all value(s)
    this.score = this.scoreService.scoreMatchFeatures(this, isolationWindow, precursorIsolationWindow, trailsWithinRTWindow);
  }

  public void matchPrecursor(ArrayList<XIC> precursors) {
    ArrayList<XIC> matchedPrecursors =
        (ArrayList<XIC>)
            subsetXICBasedOnMZ(
                precursors,
                // TODO: Configure to search different charges for now
                Utils.MassToMz(peptide.getMass(), 2));
    if (matchedPrecursors.size() > 0) {
      XIC maxIntensityPrecursor =
          matchedPrecursors.stream().max(Comparator.comparing(XIC::getMaxIntensity)).get();
      this.precursor = maxIntensityPrecursor;
      this.allRetentionTimes.add(maxIntensityPrecursor.getRtAtMaxIntensity());
    }
  }

  /*
  public void matchPrecursor(List<XIC> precursors, double precursorMZ, XIC currentPeak){
      ArrayList<XIC> matchedPrecursors = binarySearchPrecursorList(precursors, precursorMZ, currentPeak.getRtAtMaxIntensity());
      if (matchedPrecursors.size() > 0){
          XIC naivePrecursor = matchedPrecursors.get(0);
          this.precursor = new Peak(
                  naivePrecursor.getMZAtMaxIntensity(),
                  naivePrecursor.getRtAtMaxIntensity(),
                  naivePrecursor.getMaxIntensity()
          );
      }
  }*/

  public void matchIons(List<XIC> extractedPeaks, HashMap<XIC, Integer> trailToMatchedPeptidesSet) {
    // For each XIC ...
    // Assign it to one of the ArrayLists
    this.bIonsMatched = this.matchSeries(peptide.b_ions, extractedPeaks, trailToMatchedPeptidesSet);
    this.yIonsMatched = this.matchSeries(peptide.y_ions, extractedPeaks, trailToMatchedPeptidesSet);

    // Compute the +1 series
    this.bIonsPlusOneMatched = this.matchSeries(peptide.bIonPlusOne, extractedPeaks, trailToMatchedPeptidesSet);
    this.yIonsPlusOneMatched = this.matchSeries(peptide.yIonPlusOne, extractedPeaks, trailToMatchedPeptidesSet);

    // TODO: ignore the matched losses currently- they are not used
    // this.bIonH2OLossMatched = this.matchSeries(peptide.bIonH2OLoss, extractedPeaks);
    // this.bIonNH3LossMatched = this.matchSeries(peptide.bIonNH3Loss, extractedPeaks);
    // this.yIonH2OLossMatched = this.matchSeries(peptide.yIonH2OLoss, extractedPeaks);
    // this.yIonNH3LossMatched = this.matchSeries(peptide.yIonNH3Loss, extractedPeaks);
  }

  /**
   * @param mzSeries Array of ion series
   * @param possiblePeaks List of XIC organized by RETENTION TIME
   * @return
   */
  private ArrayList<XIC> matchSeries(double[] mzSeries, List<XIC> possiblePeaks, HashMap<XIC, Integer> trailToMatchedPeptidesSet) {
    ArrayList<XIC> matchedSeries = new ArrayList<>(); //pre-alloc the size of the mzSeries to avoid ensureCapacityInternal call
    matchedSeries.ensureCapacity(mzSeries.length + 1 ); // optimization fix

    for (double mz : mzSeries) {
      // TODO: Need to sort the possiblePeaks by mz again here
      // TODO: These should be sorted by m/z already, confirm this
      // possiblePeaks.sort(XIC::compareMzXIC);
      // possiblePeaks.sort(Comparator::comparingDouble(XIC::getMZAtMaxIntensity));

      List<XIC> match = subsetXICBasedOnMZ((ArrayList<XIC>) possiblePeaks, mz);
      if (!match.isEmpty()) {
        // Naively add the maximum intensity peak here
        // TODO: There may be a better match by shape here; this is just a naive computation
        XIC maxIntensityTrail = Collections.max(match,Comparator.comparingDouble(XIC::getMaxIntensity));

        // TODO: This is where you should be updating the trailToMatchedPeptidesSet
        // Want to keep track of the matched peaks, and unmatched peaks
        // Check if the key exists in the HashMap
        this.updateTrailToPeptideHashMap(maxIntensityTrail, trailToMatchedPeptidesSet);

        // build up the retention times set
        this.allRetentionTimes.add(maxIntensityTrail.getRtAtMaxIntensity());

        matchedSeries.add(maxIntensityTrail);
      } else {
        matchedSeries.add(null);
      }
    }
    return matchedSeries;
  }

  private void updateTrailToPeptideHashMap(XIC matchedTrail,HashMap<XIC, Integer> trailToMatchedPeptidesSet){
    /*
    * If the XIC is not hashed in the HashMap yet
    * Create an entry in the HashMap keyed by the XIC object
    * Set the initial Count for the trail to 1
    * ELSE
    * Update the match count by + 1
    * */
    if (!trailToMatchedPeptidesSet.containsKey(matchedTrail)){
      trailToMatchedPeptidesSet.put(matchedTrail, 1);
    } else {
      int currentCount = trailToMatchedPeptidesSet.get(matchedTrail);
      trailToMatchedPeptidesSet.replace(matchedTrail,
              currentCount + 1);
    }
  }



  private static ArrayList<XIC> binarySearchPrecursorList(
      List<XIC> precursors, double mz, double rt) {
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

    ArrayList<XIC> matchedPrecursors =
        possiblePrecursors.stream()
            .filter(
                xic -> (xic.getRtAtMaxIntensity() <= upper) && (lower <= xic.getRtAtMaxIntensity()))
            .collect(Collectors.toCollection(ArrayList::new));

    return matchedPrecursors;
  }

  // TODO: This function is incredibly unoptimized
  // TODO: 70% of all runtime comes from this subset function
  public static List<XIC> subsetXICBasedOnMZ(ArrayList<XIC> trails, double mz) {
    XIC upperMZXIC = new XIC(mz * mzUpperScalingFactor);
    XIC lowerMZXIC = new XIC(mz * mzLowerScalingFactor);

    int upperIndex = Collections.binarySearch(trails, upperMZXIC, XIC::compareMzXIC);
    int lowerIndex = Collections.binarySearch(trails, lowerMZXIC, XIC::compareMzXIC);
    upperIndex = parseInsertionPoint(upperIndex);
    lowerIndex = parseInsertionPoint(lowerIndex);

    return new ArrayList<>(trails.subList(lowerIndex, upperIndex));
  }

  private static int parseInsertionPoint(int index) {
    if (index < 0) {
      return (-1 * (index) - 1);
    } else {
      return index;
    }
  }
}
