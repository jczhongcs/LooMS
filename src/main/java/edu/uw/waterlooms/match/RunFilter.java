package edu.uw.waterlooms.match;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.Pair;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.neuralscoring.RetentionTimeMetadataTupleMap;
import edu.uw.waterlooms.service.ParameterService;
import org.json.JSONObject;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class RunFilter {

  // Peptide charge to consider
  final int FRAGMENT_ION_CHARGE = 1;
  final int PRECURSOR_CHARGE = 2;
  final String MS1_PRECURSOR_FILE =
      "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/output";

  public static ArrayList<XIC> parsePrecursorDataset(String ms1PrecursorFile, int charge)
      throws IOException {
    // TODO: Possibly Depreciated since now precursors are serialized in 0_window.ser & 0_window.tsv
    ArrayList<XIC> precursors = new ArrayList<>();

    FileInputStream inputStream = new FileInputStream(ms1PrecursorFile);
    Scanner sc = new Scanner(inputStream, "UTF-8");

    // Consume the header line
    sc.nextLine();

    // read the value(s) in
    while (sc.hasNextLine()) {
      String row = sc.nextLine();
      String[] tokens = row.split("\\s+");

      // double mz = Double.parseDouble(tokens[1]);
      // double rt = Double.parseDouble(tokens[2]);
      int z = Integer.parseInt(tokens[3]);
      if (z != charge) {
        // Skip addition if the precursor is not charge 2
        continue;
      }

      ArrayList<Double> ints = new ArrayList<>();
      ArrayList<Double> mzs = new ArrayList<>();
      ArrayList<Double> rts = new ArrayList<>();

      Double intensity_window_evg = Double.parseDouble(tokens[7]);

      XIC precursorTrail = new XIC(ints, mzs, rts);

      // Do a pseudo-normalization here to the intensity_window_evg
      precursorTrail.setMaxIntensity(precursorTrail.getMaxIntensity() / intensity_window_evg);
      precursors.add(precursorTrail);
    }

    return precursors;
  }

  public ArrayList<DbMatch> matchFromSerializedXIC(
      String xicFile,
      String serializedFileDir,
      String psmOutfile,
      String precursorFile,
      ArrayList<Peptide> candidatePeptides)
      throws IOException, ClassNotFoundException {
    // Do the matching block wise
    // Need a lookup table for the m/z range to the candidatePeptides
    // Select all peptides within the 1st window

    ArrayList<Pair> isolationWindows = new ArrayList<>();

    // Open the Precursor File
    IsolationWindow precursorWindow = null;
    FileInputStream precursorFileIn = new FileInputStream(precursorFile);
    ObjectInputStream precursorInputStream = new ObjectInputStream(precursorFileIn);
    precursorWindow = (IsolationWindow) precursorInputStream.readObject();
    precursorInputStream.close();
    precursorFileIn.close();
    ArrayList<XIC> precursors = precursorWindow.xics;
    precursors.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));




    //ArrayList<XIC> precursors = new ArrayList<>();

    // Open the xic file
    FileInputStream inputStream = new FileInputStream(xicFile);
    Scanner sc = new Scanner(inputStream, "UTF-8");

    // Read in from a slim file where the format for START-END mz is
    /*
     * START
     * START_MZ END_MZ
     * END
     * ...
     * */
    while (sc.hasNextLine()) {
      String line = sc.nextLine();
      if (line.contains("START")) {
        line = sc.nextLine();
        String[] range = line.split("\\s+");
        double mzLow = Double.parseDouble(range[0]);
        double mzHigh = Double.parseDouble(range[1]);
        Pair newWindow = new Pair(mzLow, mzHigh);
        isolationWindows.add(newWindow);
      }
    }

    ArrayList<DbMatch> peptideTrailMatches = new ArrayList<>();

    // For each isolation window...
    // Subset the candidatePeptideList
    // Read the actual IsolationWindow in
    // Match the peptides to trails
    // Output them
    int windowIndex = 1;
    for (Pair window : isolationWindows) {
      System.out.println("Matching window: " + windowIndex + "/" + isolationWindows.size());
      List<Peptide> peptideSubset =
          candidatePeptides.stream()
              .filter(
                  x -> {
                    // Filter out ALL peptides that fall between the m/z of window bound
                    double tmpMZ = Utils.MassToMz(x.mass, PRECURSOR_CHARGE);
                    double lBound = new Double(window.getL().toString());
                    double rBound = new Double(window.getR().toString());
                    return (lBound < tmpMZ && rBound > tmpMZ);
                  })
              .collect(Collectors.toList());

      // Read in the given IsolationWindow
      IsolationWindow deserializedWindowAndTrails = null;
      FileInputStream fileIn = new FileInputStream(serializedFileDir + windowIndex + "_window.ser");
      ObjectInputStream in = new ObjectInputStream(fileIn);
      deserializedWindowAndTrails = (IsolationWindow) in.readObject();
      in.close();
      fileIn.close();

      // Ensure the IsolationWindow XICs are sorted for the binary search...
      // these must be sorted by MZ
      // TODO: These xics are already sorted at L157 in Main
      deserializedWindowAndTrails.xics.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));

      for (int index = 0; index < peptideSubset.size(); index++) {

        // Output Every Thousand peptide match
        if (index % 1000 == 0) {
          System.out.println(
              "Window Index: "
                  + windowIndex
                  + ". Matching Index: "
                  + index
                  + "/"
                  + peptideSubset.size());
        }

        ArrayList<MatchedInterval> matchedIntervalsForPeptide =
            DbMatch.peptideMatch(peptideSubset.get(index), deserializedWindowAndTrails, precursorWindow, precursors);
        if (matchedIntervalsForPeptide != null && matchedIntervalsForPeptide.size() > 0) {

          // TODO: REMOVE FOR DEBUG
          // writeMatchIntervalMatchedFeaturesTSV(matchedIntervalsForPeptide, psmOutfile);

          MatchedInterval maxScoreInterval =
              Collections.max(matchedIntervalsForPeptide, Comparator.comparing(x -> x.score));
          peptideTrailMatches.add(
              new DbMatch(
                  peptideSubset.get(index).id,
                  peptideSubset.get(index).composition,
                  maxScoreInterval,
                  peptideSubset.get(index)));
        }
      }

      // TODO: This is where you should be updating the REtentionTimeMetadataTuple
      /*
       * For each element in the IsolationWindow.trailToMatchedPeptidesSet...
       * For each RetentionTime in the XIC (key)...
       * Update the RT : RetentionTimeMetadataTuple map
       * Write these values to a CSV for reading by WaterlooMS Neural Scoring
       */

      RetentionTimeMetadataTupleMap retentionTimeMetadataTupleMap = new RetentionTimeMetadataTupleMap();
      deserializedWindowAndTrails.trailToMatchedPeptidesSet.forEach(retentionTimeMetadataTupleMap::updateSetGivenTrail);
      retentionTimeMetadataTupleMap.normalizeRtTupleMapElements(deserializedWindowAndTrails.max);
      // Write the retentionTimeMetadataTupleMap to file here
      retentionTimeMetadataTupleMap.writeMapToTSV(windowIndex, psmOutfile);






      // Increment the window index
      windowIndex += 1;
    }

    // Write the result(s) to a file
    // sort result
    Collections.sort(peptideTrailMatches, Collections.reverseOrder());
    System.out.println(
        "The Number of Matched Peptides (unscored) is: " + peptideTrailMatches.size());

    writeXICResultTSV(peptideTrailMatches, psmOutfile);
    writeUnscoredPeptidesToTSV(peptideTrailMatches, psmOutfile);
    writeExperimentResultTSV(peptideTrailMatches, psmOutfile);
    return peptideTrailMatches;
  }

  /**
   * Given an ArrayList of DbMatch objects, write the XICS in JSON format to a file. Used for
   * querying results from dia_webapp.
   *
   * @param unscoredPeptides ArrayList of xic of DbMatches to be written to a file.
   * @param workingDirectory String the working directory the tsv is outputted to.
   * @throws IOException
   */
  public void writeXICResultTSV(ArrayList<DbMatch> unscoredPeptides, String workingDirectory)
      throws IOException {
    FileWriter fileWriter = new FileWriter(workingDirectory + "experiment_xic.tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
    bufferedWriter.write("index" + "\t" + "xic" + "\n");
    int idx = 1;
    for (DbMatch match : unscoredPeptides) {
      bufferedWriter.write(String.valueOf(idx) + "\t"); // index representing row_index
      bufferedWriter.write(
          jsonEncodeAllMatchedFragmentIonsBY(match.matchedRts.matchedFeatures) + "\n");
      idx += 1;
    }
    bufferedWriter.close();
    fileWriter.close();
  }

  public void writeAllMatchedIntervalsAndDummyResultTSVForAGivenPeptide(ArrayList<MatchedInterval> allMatchedIntervals, String workingDirectory, String peptide)
  throws IOException{
    FileWriter fileWriter = new FileWriter(workingDirectory + peptide + "_matchedIntervals_xic.tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
    bufferedWriter.write("index" + "\t" + "xic" + "\n");
    int idx = 1;
    for (MatchedInterval interval : allMatchedIntervals){
      bufferedWriter.write(idx + "\t");
      bufferedWriter.write(
              jsonEncodeAllMatchedFragmentIonsBY(interval.matchedFeatures) + "\n");
      idx += 1;
    }
    bufferedWriter.close();
    fileWriter.close();

    // Write the dummy result file
    fileWriter = new FileWriter(workingDirectory + peptide + "_matchedIntervals_dummyresult.tsv");
    bufferedWriter = new BufferedWriter(fileWriter);
    // this is the RT written under the peptide key
    bufferedWriter.write("index" + "\t" + "peptide" + "\t" + "RT"+ "\n");
    idx = 1;
    for (MatchedInterval interval : allMatchedIntervals){
      bufferedWriter.write(idx + "\t");
      bufferedWriter.write(peptide + "\t");
      if (interval.matchedFeatures.allRetentionTimes.size() == 0){
        bufferedWriter.write("0" + "\n");
      } else {
        bufferedWriter.write(
                interval.matchedFeatures.allRetentionTimes.get(0)
                        + "\n"); // RT TODO: Refactor to be RT at the MAX height
      }
      idx += 1;
    }
    bufferedWriter.close();
    fileWriter.close();
  }


  public void writeMatchIntervalMatchedFeaturesTSV(
      ArrayList<MatchedInterval> matchedIntervalsForPeptide, String workingDirectory)
      throws IOException {
    Peptide p = matchedIntervalsForPeptide.get(0).matchedFeatures.peptide;

    FileWriter fileWriter = new FileWriter(workingDirectory + "debug_" + p.composition + ".tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
    bufferedWriter.write(
        "rt"
            + "\t"
            + "score"
            + "\t"
            + "neutralLossScore"
            + "\t"
            + "precursorScore"
            + "\t"
            + "rtDistributionScore"
            + "\t"
            + "spScore"
            + "\t"
            + "rankScore"
            + "\t"
            + "logIntensityScore"
            + "\t"
            + "peaksPerTrailScore"
            + "\t"
            + "rtIntervalIntensityScore"
            + "\t"
            + "trailLogIntensityScore"
            + "\t"
            + "numIonsMatched"
                + "\t"
                + "rtStandardDeviation"
            + "\n");
    for (MatchedInterval interval : matchedIntervalsForPeptide) {
      MatchFeatures features = interval.matchedFeatures;
      if (features.allRetentionTimes.size() == 0) {
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\t");
        bufferedWriter.write(0 + "\n");
      } else {
        bufferedWriter.write(features.allRetentionTimes.get(0) + "\t");
        bufferedWriter.write(features.score + "\t");
        bufferedWriter.write(features.neutralLossScore + "\t");
        bufferedWriter.write(features.precursorScore + "\t");
        bufferedWriter.write(features.rtDistributionScore + "\t");
        bufferedWriter.write(features.spScore + "\t");
        bufferedWriter.write(features.rankScore + "\t");
        bufferedWriter.write(features.maxPeakLogIntensityScore + "\t");
        bufferedWriter.write(features.peaksPerTrailScore + "\t");
        bufferedWriter.write(features.rtIntervalIntensityScore + "\t");
        bufferedWriter.write(features.trailLogIntensityScore + "\t");
        bufferedWriter.write(features.numIonsMatched + "\t");
        bufferedWriter.write(String.valueOf(features.rtStandardDeviation));
        bufferedWriter.write("\n");
      }
    }
  }

  /**
   * Given an ArrayList of DbMatch objects, write the dataset to a file. Used for querying results
   * from dia_webapp.
   *
   * @param unscoredPeptides ArrayList of DbMatches to be written to a file.
   * @param workingDirectory String the working directory the tsv is outputted to.
   * @throws IOException
   */
  public void writeExperimentResultTSV(ArrayList<DbMatch> unscoredPeptides, String workingDirectory)
      throws IOException {
    FileWriter fileWriter = new FileWriter(workingDirectory + "experiment_result.tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
    bufferedWriter.write(
        "index"
            + "\t"
            + "ID"
            + "\t"
            + "peptide"
            + "\t"
            + "mz"
            + "\t"
            + "z"
            + "\t"
            + "precursor_mz"
            + "\t"
            + "ppm"
            + "\t"
            + "RT"
            + "\t"
            + "pRT"
            + "\t"
            + "iRT"
            + "\t"
            + "score"
            + "\t"
            + "quant"
            + "\t"
            + "precursor_peak_sum"
            + "\t"
            + "precursor_peak_area"
            + "\n");

    // bufferedWriter.write(jsonEncodeAllMatchedFragmentIonsBY(match.matchedRts.matchedFeatures));
    // // quantification json
    int idx = 1;
    for (DbMatch match : unscoredPeptides) {
        // Build quantification string
//      ArrayList<Double> quant = match.matchedRts.getPrecursor().getIntensities();
//      StringBuilder quantStr = new StringBuilder();
//      quantStr.append("[");
//      quant.forEach( q ->
//              quantStr.append(q).append(",")
//      );
//      quantStr.append("]");
      String quantStr = "";

      // Pass on writing the index if the matchedFeatures is null
      bufferedWriter.write(String.valueOf(idx) + "\t"); // index representing row_index
      bufferedWriter.write(match.id + "\t"); // ID
      bufferedWriter.write(match.composition + "\t"); // peptide
      bufferedWriter.write(match.getPeptide().mz + "\t"); // peptide mz at charge 1
      bufferedWriter.write(match.getPeptide().charge + "\t"); // fragment ion charge ** TODO: Refactor to precursor charge
      bufferedWriter.write(match.getPeptide().mz + "\t"); // TODO: Refactor to precusor mz
//      bufferedWriter.write(match.matchedRts.getPrecursor().getMZAtMaxIntensity() + "\t"); // precursor/observed mz
      bufferedWriter.write("10" + "\t"); // PPM TODO: Refactor
      if (match.matchedRts.getPrecursor() == null){
        bufferedWriter.write("0" + "\t");
      } else {
        bufferedWriter.write(
            match.matchedRts.getPrecursor().getRtAtMaxIntensity() + "\t"); // RT at the MAX height
      }
      bufferedWriter.write("0" + "\t"); // predicted Retention Time (pRT)
      bufferedWriter.write("0" + "\t"); // indexed Retention Time (iRT / uRT)
      bufferedWriter.write(match.matchedRts.score + "\t"); // rawScore
      bufferedWriter.write(quantStr + "\t"); // quant
      // bufferedWriter.write(jsonEncodeAllMatchedFragmentIonsBY(match.matchedRts.matchedFeatures));
      if (match.matchedRts.getPrecursor() != null) {
        bufferedWriter.write(match.matchedRts.getPrecursor().getPeakSum() + "\t"); // peak_sum of precursor
        bufferedWriter.write(match.matchedRts.getPrecursor().getPeakArea() + "\t"); // peak_area of precursor
      } else {
        bufferedWriter.write(0 + "\t"); // placeholder if not matched with a precursor
        bufferedWriter.write(0 + "\t");
      }
      bufferedWriter.write("\n");
      idx += 1;
    }
    bufferedWriter.close();
    fileWriter.close();
  }

  private ArrayList<Double> localComputationOfRetentionTimes(MatchFeatures matchFeatures){
    ArrayList<Double> allRetentionTimes = new ArrayList<>();
    int numFeatures = matchFeatures.bIonsMatched.size();
    for (int i = 0; i < numFeatures; i++) {
      XIC b = matchFeatures.bIonsMatched.get(i);
      XIC b_1 = matchFeatures.bIonsPlusOneMatched.get(i);
      XIC y = matchFeatures.yIonsMatched.get(i);
      XIC y_1 = matchFeatures.yIonsPlusOneMatched.get(i);
        if (b != null){
          allRetentionTimes.addAll(b.getRetentionTimes());
        }
        if (y != null){
          allRetentionTimes.addAll(y.getRetentionTimes());
        }
        if(b_1 != null){
          allRetentionTimes.addAll(b_1.getRetentionTimes());
        }
        if (y_1 != null){
          allRetentionTimes.addAll(y_1.getRetentionTimes());
        }
      }

    // If precursor is not null
    if (matchFeatures.precursor != null){
      allRetentionTimes.addAll(matchFeatures.precursor.getRetentionTimes());
    }

    return (ArrayList<Double>) allRetentionTimes.stream().distinct().collect(Collectors.toList());
  }


  public String jsonEncodeAllMatchedFragmentIonsBY(MatchFeatures matchedFeatures) {
    double minRT = Double.MAX_VALUE;
    double maxRT = 0;
    double maxIntensity = 0;

    ArrayList<Double> allRetentionTimes = localComputationOfRetentionTimes(matchedFeatures);

    JSONObject bIons = new JSONObject(); // RT:Intensity map
    JSONObject yIons = new JSONObject(); // RT:Intensity map
    JSONObject bMZ = new JSONObject(); // m/z map
    JSONObject yMZ = new JSONObject(); // m/z map

    JSONObject bPlusOne = new JSONObject(); // RT:Intensity map
    JSONObject yPlusOne = new JSONObject(); // RT:Intensity map
    JSONObject bPlusOneMZ = new JSONObject();
    JSONObject yPlusOneMZ = new JSONObject();

    JSONObject precursor = new JSONObject();

    JSONObject allRTs = new JSONObject();


    if (allRetentionTimes.size() > 0){
      for (int rtIdx = 0; rtIdx < allRetentionTimes.size(); rtIdx++){
        allRTs.put(String.valueOf(rtIdx), String.valueOf(allRetentionTimes.get(rtIdx)));
      }
    }

    if (matchedFeatures.precursor != null){
      JSONObject precursorTrail = new JSONObject();
      for (int precursorIdx = 0; precursorIdx < matchedFeatures.precursor.getIntensities().size(); precursorIdx ++)
      {
        double retentionTime = matchedFeatures.precursor.getRetentionTimes().get(precursorIdx);
        double intensity = matchedFeatures.precursor.getIntensities().get(precursorIdx);
        precursorTrail.put(String.valueOf(retentionTime), String.valueOf(intensity));
      }
      precursor.put(String.valueOf(Utils.MassToMz(matchedFeatures.peptide.getMass(), 2)), precursorTrail);
    }

    int numFeatures = matchedFeatures.bIonsMatched.size();
    for (int i = 0; i < numFeatures; i++) {
      bIons.put(String.valueOf(i), "");
      yIons.put(String.valueOf(i), "");
      bPlusOne.put(String.valueOf(i), "");
      yPlusOne.put(String.valueOf(i), "");


      if (matchedFeatures.bIonsPlusOneMatched.get(i) != null) {
        JSONObject feature = new JSONObject();
        for (int trailIndex = 0;
             trailIndex < matchedFeatures.bIonsPlusOneMatched.get(i).getIntensities().size();
             trailIndex++) {
          // formulate a tuple set containing rt, intensity
          double intensity = matchedFeatures.bIonsPlusOneMatched.get(i).getIntensities().get(trailIndex);
          double retentionTime =
                  matchedFeatures.bIonsPlusOneMatched.get(i).getRetentionTimes().get(trailIndex);
          if (maxIntensity < intensity) {
            maxIntensity = intensity;
          } // Update Intensity if new value is higher
          if (maxRT < retentionTime) {
            maxRT = retentionTime;
          } // Update RT if new RT is higher
          if (minRT > retentionTime) {
            minRT = retentionTime;
          } // Update RT if new RT is lower
          feature.put(String.valueOf(retentionTime), String.valueOf(intensity));
        }
        bPlusOne.put(String.valueOf(i), feature);
        bPlusOneMZ.put(String.valueOf(i), String.valueOf(matchedFeatures.bIonsPlusOneMatched.get(i).getMZAtMaxIntensity()));
      }

      if (matchedFeatures.yIonsPlusOneMatched.get(i) != null) {
        JSONObject feature = new JSONObject();
        for (int trailIndex = 0;
             trailIndex < matchedFeatures.yIonsPlusOneMatched.get(i).getIntensities().size();
             trailIndex++) {
          // formulate a tuple set containing rt, intensity
          double intensity = matchedFeatures.yIonsPlusOneMatched.get(i).getIntensities().get(trailIndex);
          double retentionTime =
                  matchedFeatures.yIonsPlusOneMatched.get(i).getRetentionTimes().get(trailIndex);
          if (maxIntensity < intensity) {
            maxIntensity = intensity;
          } // Update Intensity if new value is higher
          if (maxRT < retentionTime) {
            maxRT = retentionTime;
          } // Update RT if new RT is higher
          if (minRT > retentionTime) {
            minRT = retentionTime;
          } // Update RT if new RT is lower
          feature.put(String.valueOf(retentionTime), String.valueOf(intensity));
        }
        yPlusOne.put(String.valueOf(i), feature);
        yPlusOneMZ.put(String.valueOf(i), String.valueOf(matchedFeatures.yIonsPlusOneMatched.get(i).getMZAtMaxIntensity()));
      }

      if (matchedFeatures.bIonsMatched.get(i) != null) {
        JSONObject feature = new JSONObject();
        for (int trailIndex = 0;
            trailIndex < matchedFeatures.bIonsMatched.get(i).getIntensities().size();
            trailIndex++) {
          // formulate a tuple set containing rt, intensity
          double intensity = matchedFeatures.bIonsMatched.get(i).getIntensities().get(trailIndex);
          double retentionTime =
              matchedFeatures.bIonsMatched.get(i).getRetentionTimes().get(trailIndex);
          if (maxIntensity < intensity) {
            maxIntensity = intensity;
          } // Update Intensity if new value is higher
          if (maxRT < retentionTime) {
            maxRT = retentionTime;
          } // Update RT if new RT is higher
          if (minRT > retentionTime) {
            minRT = retentionTime;
          } // Update RT if new RT is lower
          feature.put(String.valueOf(retentionTime), String.valueOf(intensity));
        }
        bIons.put(String.valueOf(i), feature);
        bMZ.put(String.valueOf(i), String.valueOf(matchedFeatures.bIonsMatched.get(i).getMZAtMaxIntensity()));
      }
      if (matchedFeatures.yIonsMatched.get(i) != null) {
        JSONObject feature = new JSONObject();
        for (int trailIndex = 0;
            trailIndex < matchedFeatures.yIonsMatched.get(i).getIntensities().size();
            trailIndex++) {
          // formulate a tuple set containing rt, intensity
          double intensity = matchedFeatures.yIonsMatched.get(i).getIntensities().get(trailIndex);
          double retentionTime =
              matchedFeatures.yIonsMatched.get(i).getRetentionTimes().get(trailIndex);
          if (maxIntensity < intensity) {
            maxIntensity = intensity;
          } // Update Intensity if new value is higher
          if (maxRT < retentionTime) {
            maxRT = retentionTime;
          } // Update RT if new RT is higher
          if (minRT > retentionTime) {
            minRT = retentionTime;
          } // Update RT if new RT is lower
          feature.put(String.valueOf(retentionTime), String.valueOf(intensity));
        }
        yIons.put(String.valueOf(i), feature);
        yMZ.put(String.valueOf(i), String.valueOf(matchedFeatures.yIonsMatched.get(i).getMZAtMaxIntensity()));
      }
    }

    JSONObject fragmentIons = new JSONObject();
    fragmentIons.put("b", bIons);
    fragmentIons.put("y", yIons);
    fragmentIons.put("bMZ", bMZ);
    fragmentIons.put("yMZ", yMZ);

    fragmentIons.put("bPlusOne", bPlusOne);
    fragmentIons.put("yPlusOne", yPlusOne);
    fragmentIons.put("bPlusOneMZ", bPlusOneMZ);
    fragmentIons.put("yPlusOneMZ", yPlusOneMZ);

    fragmentIons.put("precursor", precursor);
    fragmentIons.put("allRTs", allRTs);
    fragmentIons.put("minRT", String.valueOf(minRT));
    fragmentIons.put("maxRT", String.valueOf(maxRT));
    fragmentIons.put("maxIntensity", String.valueOf(maxIntensity));


    return fragmentIons.toString();
  }

  /**
   * Given an ArrayList of DbMatch objects, write the score breakdown to a file.
   *
   * @param unscoredPeptides ArrayList of DbMatches to be written to a file.
   * @param workingDirectory String the working directory the tsv is outputted to.
   * @throws IOException
   */
  public void writeUnscoredPeptidesToTSV(
      ArrayList<DbMatch> unscoredPeptides, String workingDirectory) throws IOException {
    FileWriter fileWriter = new FileWriter(workingDirectory + "experiment_scores.tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
    // Write the Header
    bufferedWriter.write(
        "index"
            + "\t"
            + "id"
            + "\t"
            + "sequence"
            + "\t"
            + "rawScore"
            + "\t"
            + "fragmentIonCharge"
            + "\t"
            + "logIntensityScore"
            + "\t"
            + "spScore"
            + "\t"
            + "neutralLossScore"
            + "\t"
            + "rtDistributionScore"
            + "\t"
            + "rankScore"
            + "\t"
            + "precursorScore"
            + "\t"
            + "peaksPerTrailScore"
            + "\t"
            + "trailLogIntensityScore"
            + "\t"
            + "rtIntervalIntensityScore"
            + "\t"
            + "numIonsMatched"
                + "\t"
                + "rtStandardDeviation"
            + "\n");

    int idx = 1;
    for (DbMatch match : unscoredPeptides) {
      // Output the information to a file
      bufferedWriter.write(String.valueOf(idx) + "\t"); // index representing row_index
      bufferedWriter.write(match.id + "\t");
      bufferedWriter.write(match.composition + "\t");
      bufferedWriter.write(Double.toString(match.matchedRts.score) + "\t");
      bufferedWriter.write(Integer.toString(this.FRAGMENT_ION_CHARGE) + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.maxPeakLogIntensityScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.spScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.neutralLossScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.rtDistributionScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.rankScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.precursorScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.peaksPerTrailScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.trailLogIntensityScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.rtIntervalIntensityScore + "\t");
      bufferedWriter.write(match.matchedRts.matchedFeatures.numIonsMatched + "\t");
      bufferedWriter.write(String.valueOf(match.matchedRts.matchedFeatures.rtStandardDeviation));
      bufferedWriter.write("\n");
      idx += 1;
    }
    bufferedWriter.close();
    fileWriter.close();
  }

  /**
   * Writes an ArrayList of Peptide objects to a TSV file, with the b and y ions separated by
   * commas.
   *
   * @param candidatePeptides ArrayList<Peptide> of peptides to be written.
   * @param workingDirectory String of directory the file is being written to.
   * @throws IOException
   */
  public void writePeptidesToTSV(ArrayList<Peptide> candidatePeptides, String workingDirectory)
      throws IOException {
    FileWriter fileWriter = new FileWriter(workingDirectory + "candidatePeptides.tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
    // Write the header
    bufferedWriter.write(
        "id"
            + "\t"
            + "sequence"
            + "\t"
            + "mass"
            + "\t"
            + "mz"
            + "\t"
            + "charge"
            + "\t"
            + "y_ions"
            + "\t"
            + "b_ions"
            + "\n");

    for (Peptide p : candidatePeptides) {
      bufferedWriter.write(p.id + "\t");
      bufferedWriter.write(p.composition + "\t");
      bufferedWriter.write(p.mass + "\t");
      bufferedWriter.write(p.mz + "\t");
      bufferedWriter.write(p.charge + "\t");

      // Write y_ions
      StringBuilder stringBuilder = new StringBuilder();
      for (Double yion : p.y_ions) {
        stringBuilder.append(yion).append(",");
      }
      bufferedWriter.write(stringBuilder.toString() + "\t");

      // Write b_ions
      stringBuilder = new StringBuilder();
      for (Double bion : p.b_ions) {
        stringBuilder.append(bion).append(",");
      }
      bufferedWriter.write(stringBuilder.toString() + "\n");
    }

    bufferedWriter.close();
    fileWriter.close();
  }
}
