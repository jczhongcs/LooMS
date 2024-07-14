package edu.uw.waterlooms.utility;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.Pair;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.match.*;
import org.junit.Test;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.*;
import java.util.stream.Collectors;

public class GeneratePeptideTrailIntegrationTest {
  String r0 = "1";
  String q0 = "1";
  /*String R01_Q1_INPUT =
      "/home/jia/Documents/Code/waterlooms/wms_neural_scoring/data/development_data/r0"
          + r0
          + "_q0"
          + q0
          + "_mzid.tsv";*/
  String R01_Q1_INPUT = "/Users/jianzhong/Documents/uwaterloo/dia_data_reading-looms-refactor/dia_data_reading/data/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_Q1_stripped_for_training_mzid.tsv";


//  String SNE_INPUT =
//      "/home/jia/Documents/Code/waterlooms/wms_neural_scoring/data/development_data/sne_stripped.tsv";
//  String MUS_INPUT =
//      "/home/jia/Documents/Code/waterlooms/wms_neural_scoring/data/development_data/mus_stripped.tsv";
  final String XIC_OUTPUT_DIR = "/Users/jianzhong/Documents/uwaterloo/dia_data_reading-looms-refactor/dia_data_reading/data/output/";
  final String SERIALIZED_DATA_DIR =
      "/Users/jianzhong/Documents/uwaterloo/dia_data_reading-looms-refactor/dia_data_reading/data/output/";

  private HashMap<String, Double> parseSNESubset(String sneSubsetFile) {
    HashMap<String, Double> peptideAndRT = new HashMap<>();

    try {
      FileInputStream inputStream = new FileInputStream(sneSubsetFile);
      Scanner sc = new Scanner(inputStream, "UTF-8");

      // Consume the header
      sc.nextLine();
      while (sc.hasNextLine()) {
        String line = sc.nextLine();
        String[] tokens = line.split("\\s+");
        peptideAndRT.put(tokens[0], Double.valueOf(tokens[1]));
      }

      sc.close();
      inputStream.close();
    } catch (Exception $e) {
      return peptideAndRT;
    }
    return peptideAndRT;
  }

  private HashMap<String, Double> parseMSGFSubset(String msgfSubsetFile) {
    HashMap<String, Double> peptideAndRT = new HashMap<>();

    try {
      FileInputStream inputStream = new FileInputStream(msgfSubsetFile);
      Scanner sc = new Scanner(inputStream, "UTF-8");

      // Consume the header
      sc.nextLine();
      while (sc.hasNextLine()) {
        String line = sc.nextLine();
        String[] tokens = line.split("\\s+");
        peptideAndRT.put(tokens[19], Double.valueOf(tokens[18]));
      }

      sc.close();
      inputStream.close();
    } catch (Exception $e) {
      return peptideAndRT;
    }
    return peptideAndRT;
  }

  // TODO: This is used to generate the XIC data for Xiangyuan's quantification component
//  @Test
  public void generateXICsUsingMSGF_RTs() throws IOException {
    // Given
    // Uncomment for MSGF+ Results
    HashMap<String, Double> r01 = parseMSGFSubset(this.R01_Q1_INPUT);
    // Uncomment for SNE Results
    //HashMap<String, Double> r01 = parseSNESubset(this.SNE_INPUT);

    ArrayList<Peptide> candidatePeptides = new ArrayList<>();

    // For each PEPTIDE_SEQUENCE:RT element ...
    r01.forEach(
        (sequence, rt) -> {
          // Generate the Fragment Ions
          Peptide peptide = new Peptide("1%fdr", (String) sequence);
          // Append them to the list
          candidatePeptides.add(peptide);
        });

    ArrayList<DbMatch> allMatches = new ArrayList<>();

    ArrayList<Pair> IsolationWindowRanges =
        this.parseWindowRanges(
            "/Users/jianzhong/Documents/uwaterloo/dia_data_reading-looms-refactor/dia_data_reading/data/output/isolationWindowRanges.out");

    /*
    ArrayList<XIC> precursors =
    RunFilter.parsePrecursorDataset(
        "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/r0"
            + r0
            + "_dia_data.mzXML.precursors",
        2); */

//    /*
    ArrayList<XIC> precursors = new ArrayList<>();
    // Open the Precursor File
    IsolationWindow precursorWindow = null;
    try {
      FileInputStream precursorFileIn =
          new FileInputStream(
              "/Users/jianzhong/Documents/uwaterloo/dia_data_reading-looms-refactor/dia_data_reading/data/output/0_window.ser");
      ObjectInputStream precursorInputStream = new ObjectInputStream(precursorFileIn);
      precursorWindow = (IsolationWindow) precursorInputStream.readObject();
      precursorInputStream.close();
      precursorFileIn.close();
      precursors = precursorWindow.xics;
      precursors.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));
    } catch (ClassNotFoundException $e) {
    }
//     */

    RunFilter runFilter = new RunFilter();
    int windowIndex = 1;
    for (Pair window : IsolationWindowRanges) {
      System.out.println("Matching window: " + windowIndex + "/" + IsolationWindowRanges.size());
      List<Peptide> peptideSubset =
          candidatePeptides.stream()
              .filter(
                  x -> {
                    // Filter out ALL peptides that fall between the m/z of window bound
                    double tmpMZ = Utils.MassToMz(x.mass, 2);
                    double lBound = new Double(window.getL().toString());
                    double rBound = new Double(window.getR().toString()) - 1;
                    return (lBound < tmpMZ && rBound >= tmpMZ);
                  })
              .collect(Collectors.toList());

      System.out.println("There are " + peptideSubset.size() + " peptides in this window.");
      if (peptideSubset.size() == 0) {
        windowIndex += 1;
        continue;
      }
      // Read in the given IsolationWindow
      IsolationWindow deserializedWindowAndTrails = parseSerializedIsolationWindowFile(windowIndex);
      assert deserializedWindowAndTrails != null;
      deserializedWindowAndTrails.xics.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));
      for (int index = 0; index < peptideSubset.size(); index++) {
        ArrayList<MatchedInterval> matchedIntervalsForPeptide =
            DbMatch.peptideMatch(peptideSubset.get(index), deserializedWindowAndTrails, precursorWindow, precursors);
        if (matchedIntervalsForPeptide.size() > 0) {
          int finalIndex = index;
          MatchedInterval closestRTInterval =
              Collections.min(
                  matchedIntervalsForPeptide,
                  Comparator.comparingDouble(
                      x -> {
                        double initialRT = 0;
                        if (x.getMatchedFeatures().allRetentionTimes.size() == 0) {
                          initialRT = 0;
                        } else {
                          initialRT = x.getMatchedFeatures().allRetentionTimes.get(0);
                        }
                        double targetRT = r01.get(peptideSubset.get(finalIndex).composition);
                        return Math.abs(initialRT - targetRT);
                      }));
          allMatches.add(
              new DbMatch(
                  peptideSubset.get(index).id,
                  peptideSubset.get(index).composition,
                  closestRTInterval,
                  peptideSubset.get(index)));
        }
      }
      // TODO: REMOVE THIS
      // runFilter.writeXICResultTSV(allMatches, this.XIC_OUTPUT_DIR + "r01_q01_xic_");

      // System.exit(0);
      windowIndex += 1;
    }

    try {
      //* Uncomment for Writing for MSGF+
      runFilter.writeXICResultTSV(
          allMatches, this.XIC_OUTPUT_DIR + "r0" + r0 + "_q0" + q0 + "_xic_");
      runFilter.writeUnscoredPeptidesToTSV(
          allMatches, this.XIC_OUTPUT_DIR + "r0" + r0 + "_q0" + q0 + "_xic_");
      runFilter.writeExperimentResultTSV(
          allMatches, this.XIC_OUTPUT_DIR + "r0" + r0 + "_q0" + q0 + "_xic_");
      //*/
      /* Uncomment for Writing for SNE
      runFilter.writeXICResultTSV(allMatches, this.XIC_OUTPUT_DIR + "r01_sne_xic_");
      runFilter.writeUnscoredPeptidesToTSV(allMatches, this.XIC_OUTPUT_DIR + "r01_sne_xic_");
      runFilter.writeExperimentResultTSV(allMatches, this.XIC_OUTPUT_DIR + "r01_sne_xic_");
      // */
    } catch (Exception ignored) {
    }
    // When
    // Then
  }
  /**
   * Given a integer representing the window to be searched, deserialize the IsolationWindow.
   *
   * @param windowIndex int representing the window to be searched.
   * @return an unserialized IsolationWindow object.
   */
  private IsolationWindow parseSerializedIsolationWindowFile(int windowIndex) {
    IsolationWindow deserializedIsolationWindowTrails = null;

    try {
      FileInputStream fileIn =
          new FileInputStream(this.SERIALIZED_DATA_DIR + windowIndex + "_window.ser");
      ObjectInputStream in = new ObjectInputStream(fileIn);
      deserializedIsolationWindowTrails = (IsolationWindow) in.readObject();
      in.close();
      fileIn.close();

    } catch (IOException | ClassNotFoundException $e) {
      return null;
    }
    return deserializedIsolationWindowTrails;
  }

  private ArrayList<Pair> parseWindowRanges(String file) throws FileNotFoundException {
    try {
      ArrayList<Pair> isolationWindows = new ArrayList<>();
      FileInputStream inputStream = new FileInputStream(file);
      Scanner sc = new Scanner(inputStream, "UTF-8");
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
      return isolationWindows;
    } catch (FileNotFoundException $e) {
      System.out.println("error");
    }
    return new ArrayList<>();
  }
}
