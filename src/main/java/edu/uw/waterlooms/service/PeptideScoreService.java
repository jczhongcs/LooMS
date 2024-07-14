package edu.uw.waterlooms.service;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.entity.XICEntry;
import edu.uw.waterlooms.match.MatchFeatures;
import edu.uw.waterlooms.match.Peak;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;
import java.util.stream.Collectors;

import static java.lang.Math.log10;
import static java.lang.Math.max;

public class PeptideScoreService {
  private String ms1FeatureFile;
  private ArrayList<XICEntry> features;

  private static final Double RT_TOLERANCE_MINUTES = 0.16;

  /** Empty Base Constructor. */
  public PeptideScoreService() {}

  /**
   * PeptideScoreService Constructor
   *
   * @param ms1FeatureFile String representing the filtered MS1 feature file.
   */
  public PeptideScoreService(String ms1FeatureFile) {
    this.ms1FeatureFile = ms1FeatureFile;
    this.features = new ArrayList<>();
  }

  /** @return ArrayList<XICEntry> of MS1 features. */
  public ArrayList<XICEntry> getFeatures() {
    return this.features;
  }

  /** Helper function to populate the ArrayList of XICEntries from the TSV */
  public void populateMS1Features() throws IOException {
    /*
    Open TSV
    Read line by line
    Add to the MS1 feature file.
    */
    FileInputStream inputStream = new FileInputStream(this.ms1FeatureFile);
    Scanner scanner = new Scanner(inputStream, "UTF-8");

    // Iterate through the lines to generate the XICEntries
    scanner.nextLine(); // Consume the header

    while (scanner.hasNextLine()) {
      String line = scanner.nextLine();
      // Delimit line by TAB character and populate the XICEntry
      String[] tokens = line.split("\\t");
      XICEntry xic =
          new XICEntry(
              Double.parseDouble(tokens[1]), // mz
              Double.parseDouble(tokens[2]), // rt
              Double.parseDouble(tokens[0])); // intensity TODO: DEPRECIATED CURRENTLY
      // Append this to the feature list
      // Should already be sorted by retentionTime
      this.features.add(xic);
    }

    // Finally close the stream(s)
    scanner.close();
    inputStream.close();
  }

  /**
   * Compute the score contribution from CS482.
   *
   * @param x a normalized value, typically (intensity / maxIntensity)
   * @return double score
   */
  public static double calculateScoreContribution(double x) {
    // Implement score f'n from class
    double scoreContribution = 0d;
    if (x > 0.01) {
      scoreContribution = log10(100 * x);
    }
    return max(0d, scoreContribution);
  }

  public double scoreMatchFeatures(
      MatchFeatures matchFeatures,
      IsolationWindow isolationWindow,
      IsolationWindow precursorIsolationWindow,
      List<XIC> trailsWithinRTWindow) {
    boolean consecutive = false;
    double reward = 0;
    int m = 0;

    double spScore = 0;
    double rankScore = 0;
    double neutralLossScore = 0;
    double precursorScore = 0;
    double rtDistributionScore = 0;
    double peakLogIntensityScore = 0;
    double peaksPerTrailScore = 0;
    double trailLogIntensityScore = 0;
    double rtStandardDeviationScore = 0;

    // peaksPerTrailScore is just spScore adjusted by the number of trails that contribute to that
    // peak
    // theory is that the trails with more peaks should be weighted more heavily

    // Iterate across matchFeatures
    for (int i = 0; i < matchFeatures.bIonsMatched.size(); i++) {
      XIC bIon = matchFeatures.bIonsMatched.get(i);
      XIC yIon = matchFeatures.yIonsMatched.get(i);

      trailLogIntensityScore += calculateTrailLogIntensityScoreContribution(bIon, isolationWindow);
      trailLogIntensityScore += calculateTrailLogIntensityScoreContribution(yIon, isolationWindow);
      double bTmp = calculateScoreContribution(normalizeAgainstMaxIntensity(bIon, isolationWindow));
      double yTmp = calculateScoreContribution(normalizeAgainstMaxIntensity(yIon, isolationWindow));
      peakLogIntensityScore += bTmp;
      peakLogIntensityScore += yTmp;
    }
    precursorScore += calculatePrecursorContribution(matchFeatures.precursor, precursorIsolationWindow);

    double cumulativeScore = trailLogIntensityScore + peakLogIntensityScore + precursorScore;

    matchFeatures.rankScore = 0;
    matchFeatures.neutralLossScore = 0;
    matchFeatures.precursorScore = 0;
    matchFeatures.rtDistributionScore = 0;
    matchFeatures.spScore = 0;
    matchFeatures.maxPeakLogIntensityScore = peakLogIntensityScore;
    matchFeatures.peaksPerTrailScore = 0;
    matchFeatures.rtIntervalIntensityScore = 0;
    matchFeatures.trailLogIntensityScore = trailLogIntensityScore;
    matchFeatures.rtStandardDeviation = 0;
    return cumulativeScore;

    /*
      XIC bIonH2OLoss = matchFeatures.bIonH2OLossMatched.get(i);
      XIC bIonNH3Loss = matchFeatures.bIonNH3LossMatched.get(i);
      XIC yIonH2OLoss = matchFeatures.yIonH2OLossMatched.get(i);
      XIC yIonNH3Loss = matchFeatures.yIonNH3LossMatched.get(i);

      if (bIon != null){
        peaksPerTrailScore += (bIon.getRetentionTimes().size() / 5d); // TODO: change 5d to a parameter
      }
      if (yIon != null){
        peaksPerTrailScore +=  (yIon.getRetentionTimes().size() / 5d);
      }
      spScore += bTmp;
      spScore += yTmp;
      //spScore += calculateScoreContribution(normalizeAgainstMaxIntensity(yIon, isolationWindow));
      //spScore += calculateScoreContribution(normalizeAgainstMaxIntensity(bIonH2OLoss, isolationWindow));
      //spScore += calculateScoreContribution(normalizeAgainstMaxIntensity(yIonH2OLoss, isolationWindow));
      //spScore += calculateScoreContribution(normalizeAgainstMaxIntensity(bIonNH3Loss, isolationWindow));
      //spScore += calculateScoreContribution(normalizeAgainstMaxIntensity(yIonNH3Loss, isolationWindow));
      // Add the reward for consecutive ion series matching
      if (bIon != null || yIon != null){
        consecutive = true;
        m += 1;
      } else {
        consecutive = false;
      }
      if (consecutive){
        reward += 0.075;
      }

      // Compute RANK Score
      rankScore += calculateRankContribution(bIon, isolationWindow);
      rankScore += calculateRankContribution(yIon, isolationWindow);
      //rankScore += calculateRankContribution(bIonH2OLoss, isolationWindow);
      //rankScore += calculateRankContribution(yIonH2OLoss, isolationWindow);
      //rankScore += calculateRankContribution(bIonNH3Loss, isolationWindow);
      //rankScore += calculateRankContribution(yIonNH3Loss, isolationWindow);

      // Compute Neutral Loss Score
      // TODO: Neutral Loss only works if you have a matching bIon
      neutralLossScore += calculateNeutralLossContribution(bIon, bIonH2OLoss, bIonNH3Loss, isolationWindow);
      neutralLossScore += calculateNeutralLossContribution(yIon, yIonH2OLoss, yIonNH3Loss, isolationWindow);
    }
    // Compute Precursor Score

    // Compute RT Distribution Score
    rtDistributionScore += calculateRTDistributionContribution(matchFeatures);

    // Calculate the final Intensity Score (Sp) score in Comet
    spScore = (spScore * m * (1 + reward)) / (matchFeatures.bIonsMatched.size() - 1);

    rtStandardDeviationScore = calculateSD(matchFeatures.allRetentionTimes);

    // Sum up all the score(s)
    // This allows you to WEIGHT the contribution(s) differently)
    // WRITE THIS TO A FILE FOR SAVING
    // double cumulativeScore = spScore + rankScore + neutralLossScore + precursorScore + rtDistributionScore + peaksPerTrailScore;
    double cumulativeScore = spScore + trailLogIntensityScore + rankScore + neutralLossScore + precursorScore + peaksPerTrailScore - rtDistributionScore - rtStandardDeviationScore;

    matchFeatures.rankScore = rankScore;
    matchFeatures.neutralLossScore = neutralLossScore;
    matchFeatures.precursorScore = precursorScore;
    matchFeatures.rtDistributionScore = rtDistributionScore;
    matchFeatures.spScore = spScore;
    matchFeatures.maxPeakLogIntensityScore = peakLogIntensityScore;
    matchFeatures.peaksPerTrailScore = peaksPerTrailScore;
    matchFeatures.rtIntervalIntensityScore = calculateRTIntervalIntensityContribution(matchFeatures, trailsWithinRTWindow);
    matchFeatures.trailLogIntensityScore = trailLogIntensityScore;
    matchFeatures.rtStandardDeviation = rtStandardDeviationScore;
    return cumulativeScore;
    */
  }

  public static double calculateRTIntervalIntensityContribution(
      MatchFeatures matchFeatures, List<XIC> trailsWithinRTWindow) {
    /*
    Extract all XICs within the RT interval
    Extract normalize matched RTs to the intensity
    */

    double rSum = 0d;
    for (XIC trail : trailsWithinRTWindow) {
      rSum += trail.getIntensities().stream().mapToDouble(Double::doubleValue).sum();
    }

    double intensitySum = 0d;
    // TODO: optimize here
    for (int i = 0; i < matchFeatures.bIonsMatched.size(); i++) {
      XIC bIons = matchFeatures.bIonsMatched.get(i);
      XIC yIons = matchFeatures.yIonsMatched.get(i);
      if (bIons != null) {
        intensitySum += bIons.getIntensities().stream().mapToDouble(Double::doubleValue).sum();
        matchFeatures.numIonsMatched += 1;
      }
      if (yIons != null) {
        intensitySum += yIons.getIntensities().stream().mapToDouble(Double::doubleValue).sum();
        matchFeatures.numIonsMatched += 1;
      }
    }

    if (rSum != 0) {
      return intensitySum / rSum;
    } else {
      return 0;
    }
  }

  public static double calculateRTDistributionContribution(MatchFeatures matchFeatures) {
    if (matchFeatures.allRetentionTimes.size() == 0) {
      return 0;
    }
    double maxRT = Collections.max(matchFeatures.allRetentionTimes);
    double minRT = Collections.min(matchFeatures.allRetentionTimes);

    return (maxRT - minRT);
  }

  public static double calculatePrecursorContribution(
      XIC precursor, IsolationWindow isolationWindow) {
    if (precursor == null) {
      return 0;
    }

    // Precursor is already normalized
    // TODO: FIX THIS
    return calculateScoreContribution(normalizeAgainstMaxIntensity(precursor, isolationWindow));
  }

  public static double calculateNeutralLossContribution(
      XIC baseTrail, XIC H2OTrail, XIC NH3Trail, IsolationWindow isolationWindow) {
    double score = 0;
    // Contribute to the score IFF base Peak exists
    if (baseTrail == null) {
      return 0;
    }
    if (H2OTrail != null) {
      // Calculate here
      score += calculateScoreContribution(normalizeAgainstMaxIntensity(H2OTrail, isolationWindow));
    }
    if (NH3Trail != null) {
      score += calculateScoreContribution(normalizeAgainstMaxIntensity(NH3Trail, isolationWindow));
    }
    return score;
  }

  public static double calculateRankContribution(XIC trail, IsolationWindow isolationWindow)
        // TODO: This consumes about 14.5% of the runtime of the scoring, perhaps need to do this
        // binary search in a better way
      {
    if (trail == null) {
      return 0;
    }
    double rank =
        Collections.binarySearch(
            isolationWindow.rtToIntensitySet.get(trail.getRtAtMaxIntensity()),
            trail.getMaxIntensity());
    return rank / isolationWindow.rtToIntensitySet.get(trail.getRtAtMaxIntensity()).size();
  }

  public static double calculateSD(ArrayList<Double> doubleArray) {
    Double sum = doubleArray.stream().mapToDouble(Double::doubleValue).sum();
    Double mean = sum / doubleArray.size();
    Double sd = 0d;

    for (Double num : doubleArray) {
      sd += Math.pow(num - mean, 2);
    }

    return Math.sqrt(sd / doubleArray.size());
  }

  public static double calculateTrailLogIntensityScoreContribution(
      XIC trail, IsolationWindow isolationWindow) {
    if (trail == null) {
      return 0;
    }
    // for each trail peak
    // get the intensity
    // get the rt
    // normalize to max intensity
    // take summation of these
    // calculate score contribution
    double trailScore = 0d;

    // Iterate across the set of trail(s)
    for (int idx = 0; idx < trail.getIntensities().size(); idx++) {
      trailScore +=
          trail.getIntensities().get(idx)
              / isolationWindow.max.get(trail.getRetentionTimes().get(idx));
    }

    return calculateScoreContribution(trailScore);
  }

  private static double normalizeAgainstMaxIntensity(XIC trail, IsolationWindow isolationWindow) {
    if (trail == null) {
      return 0;
    }
    double rt = trail.getRtAtMaxIntensity();
    // System.out.println("rt: " + rt);
    double intensity = trail.getMaxIntensity();
    // System.out.println("intensity: " + intensity);
    // System.out.println("isolationWindow.max.get(rt): " + isolationWindow.max.get(rt));
    double intensityNormalization = isolationWindow.max.get(rt);

    return (intensity / intensityNormalization);
  }
}
