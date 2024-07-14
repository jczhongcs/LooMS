package edu.uw.waterlooms.match;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.Pair;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.ms1.FeatureDetect;
import edu.uw.waterlooms.msutil.OpenMzxml;
import edu.uw.waterlooms.neuralscoring.RetentionTimeMetadataTupleMap;
import org.junit.Test;

import java.io.*;
import java.nio.file.FileSystems;
import java.util.*;
import java.util.stream.Collectors;

public class RunFilterIntegrationTest {

  // TODO: Change this to the resource(s) folder if possible
  final String WORKING_DIRECTORY = "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/debug_output/";
  final String MS1_PRECURSOR_FILE =
      "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/r01_dia_data.precursors";

  final String XIC_WINDOW_PAIRFILE =
      "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/output/isolationWindowRanges.out";
  final String FASTA_FILE = //"/home/jia/Documents/Code/waterlooms/dia_data_reading/data/database_withdecoy.fasta";
      //"/home/jia/Documents/Code/waterlooms/dia_data_reading/data/top_10_with_decoys.fasta";
      // "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/top_100_with_decoys.fasta";
      "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/top_1000_with_decoys.fasta";
  final int MINIMUM_PEPTIDE_LENGTH = 7;

  final String PEPTIDE_CANDIDATE = "SILFVPTSAPR";

  final String SERIALIZED_DATA_DIR =
      "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/output/";
  final int FRAGMENT_ION_CHARGE = 1;

  /**
   * Reads a slim .out file returning the Upper and Lower bound mz value(s) for all Isolation
   * Windows.
   *
   * @return ArrayList<Pair> of the Upper and Lower bound(s)
   */
  private ArrayList<Pair> parseIsolationWindows() {
    ArrayList<Pair> isolationWindows = new ArrayList<>();

    try {
      FileInputStream inputStream = new FileInputStream(this.XIC_WINDOW_PAIRFILE);
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
    } catch (IOException $e) {
      return isolationWindows; // Return empty copy
    }
    return isolationWindows;
  }

  /**
   * Given a mz and the set of IsolationWindows, return the window the mz falls into.
   *
   * @param windows ArrayList of IsolationWindows to be searched.
   * @param mz double MZ to be searched
   * @return int representing the window
   */
  private int findWindowGivenMZ(ArrayList<Pair> windows, double mz) {
    // The Serialized objects are indexed from 1
    int windowIndex = 1;

    for (int index = 0; index < windows.size(); index++) {
      double lower = new Double(windows.get(index).getL().toString());
      double upper = new Double(windows.get(index).getR().toString());

      if (lower <= mz && mz <= upper) {
        return windowIndex;
      }
      windowIndex += 1;
    }

    return windowIndex;
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

  // Integration test that is primarily for debug purposes.
  // Will fail if run in production since the serialized windows do not exist.
//  @Test
  public void matchACandidateSubset() throws IOException {
    ArrayList<Pair> isolationWindows = parseIsolationWindows();
    ArrayList<Peptide> topScoringCandidates = new ArrayList<>();

    // Build up the Peptide to be matched
    // TODO: These are candidates

    /*
    Peptide peptide = new Peptide("11_target", "SILFVPTSAPR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "APDFVFYAPR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "SSEHINEGETAMLVCK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "EAAENSLVAYK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "VKPYLPQICGTVLWR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "MTDQEAIQDLWQWR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "AEAGDNLGALVR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "GPLQSVQVFGR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "TDLEKDIISDTSGDFR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_target", "HSQFIGYPITLFVEK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    */

    String candidate = "DLYANTVLSGGTTMYPGIADR";
    ArrayList<Peptide> singlePeptide = new ArrayList<>();
    Peptide singleCandidate = new Peptide("first_common_target", candidate);
    singlePeptide.add(singleCandidate);
    singleCandidate.GenerateIons();

    /*
    edu.ucsd.msjava.msutil.Peptide ucsdPeptide = AminoAcidSet.getStandardAminoAcidSet().getPeptide(candidate);

    //ArrayList<Composition> bIons = ucsdPeptide.toCumulativeCompositionSequence(true, new Composition(0, 1, 0, 0, 0));
    double mass = ucsdPeptide.getMass();
    double accurateMass = ucsdPeptide.getAccurateMass();
    double nominalMass = ucsdPeptide.getNominalMass();
    double parentMass = ucsdPeptide.getParentMass();

    int tmpX = 0;*/



    /*
    Peptide peptide = new Peptide("11_decoy", "EEDEGESDNR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "DPCGSSPSCGR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "SPGCDEGQDSK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "SEDHCSESTK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "DGEDEQEETK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "GDPDCSASAGSR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "YWCSFTCR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "GPSPCDGTTCK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "TEDFHGACSR");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    peptide = new Peptide("11_decoy", "SHEGDTEYNK");
    peptide.GenerateIons();
    topScoringCandidates.add(peptide);
    */







    // Select the IsolationWindow the precursor falls within at z = 2
    Double mzAtChargeTwo = Utils.MassToMz(singleCandidate.mass, 2);

    int windowIndex = findWindowGivenMZ(isolationWindows, mzAtChargeTwo);
    IsolationWindow isolationWindow = parseSerializedIsolationWindowFile(windowIndex);

    // Need to sort IsolationWindow by mz first
    assert isolationWindow != null;
    isolationWindow.xics.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));

    final String MS1_PRECURSOR_FILE =
        "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/output/precursors.tsv";
    final int PRECURSOR_CHARGE = 2;

    // when window Index == 11, there was an interesting bug (not bug, if low enough score it's
    // negative
    //ArrayList<XIC> precursors = parsePrecursorDataset(MS1_PRECURSOR_FILE, PRECURSOR_CHARGE);
    //precursors.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));
    ArrayList<XIC> precursors = new ArrayList<>();

    ArrayList<DbMatch> allMatches = new ArrayList<>();
    for (Peptide p : singlePeptide){
      String MATCHED_INTERVAL_OUTPUT_DIR =
              "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/debug_output/"
                      + p.composition
                      + "_matchIntervals.tsv";
      ArrayList<MatchedInterval> matchedIntervalsForPeptide =
              DbMatch.peptideMatch(p, isolationWindow, isolationWindow, precursors); // TODO: This is a bug, the 2nd isolationWindow should be PRECURSOR isolation window
      assert matchedIntervalsForPeptide != null;
      MatchedInterval maxScoreInterval = Collections.max(matchedIntervalsForPeptide, Comparator.comparing(x -> x.score));
      allMatches.add(
              new DbMatch(
                      p.id,
                      p.composition,
                      maxScoreInterval,
                      p)
      );

      DbMatch.writeMatchWindows(matchedIntervalsForPeptide, MATCHED_INTERVAL_OUTPUT_DIR);
    }
  }


  /**
   * Read in the pre-determined FASTA file.
   *
   * @return ArrayList<Peptide> of CandidatePeptides.
   */
  private ArrayList<Peptide> parseFastaFile(int minimumPeptideLength) {
    ArrayList<Genome> allPeptides = null;
    ArrayList<Peptide> candidatePeptides = null;
    try {
      allPeptides = FastaFile.ReadFile(this.FASTA_FILE);
      candidatePeptides =
          allPeptides.stream()
              .map(x -> x.composition)
              .flatMap(Collection::stream)
              .distinct()
              .sorted(Comparator.comparingDouble(Peptide::getMz))
              .collect(Collectors.toCollection(ArrayList::new));
      candidatePeptides.removeIf(x -> x.composition.length() < minimumPeptideLength);
    } catch (IOException $e) {
      return null;
    }

    return candidatePeptides;
  }

  /**
   * Given the MS1 Precursors computed by Xiangyuan's component
   *
   * @param ms1PrecursorFile
   * @return an ArrayList of Precursor(s) with precursor charge filtering.
   * @throws IOException
   */
  public ArrayList<XIC> parsePrecursorDataset(String ms1PrecursorFile, int charge)
      throws IOException {
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

      // These are the isotopic distributions...
      mzs.add(Double.parseDouble(tokens[15]));
      mzs.add(Double.parseDouble(tokens[16]));
      mzs.add(Double.parseDouble(tokens[17]));
      mzs.add(Double.parseDouble(tokens[18]));

      ints.add(Double.parseDouble(tokens[19]));
      ints.add(Double.parseDouble(tokens[20]));
      ints.add(Double.parseDouble(tokens[21]));
      ints.add(Double.parseDouble(tokens[22]));

      rts.add(Double.parseDouble(tokens[23]));
      rts.add(Double.parseDouble(tokens[24]));
      rts.add(Double.parseDouble(tokens[25]));
      rts.add(Double.parseDouble(tokens[26]));

      Double intensity_window_evg = Double.parseDouble(tokens[7]);

      XIC precursorTrail = new XIC(ints, mzs, rts);

      // Do a pseudo-normalization here to the intensity_window_evg
      precursorTrail.setMaxIntensity(precursorTrail.getMaxIntensity() / intensity_window_evg);
      precursors.add(precursorTrail);
    }

    return precursors;
  }

//  @Test
  public void matchAFullIsolationWindow() throws IOException {
    // Takes about 3 minutes for window 11
    // This test assumes that Xiangyuan's MS1 feature detection and selection has been performed
    // Given a full set of CandidatePeptide(s)
    // Lets say window 11

    final int PRECURSOR_CHARGE = 2;
    ArrayList<XIC> precursors = new ArrayList<>();
    //ArrayList<XIC> precursors = parsePrecursorDataset(this.MS1_PRECURSOR_FILE, PRECURSOR_CHARGE);
    // Sort the Precursor set by MZ for binary searching
    //precursors.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));
    // precursors = new ArrayList(); // TODO: Uncomment this for NO precursor scoring

    int windowIndex = 11;
    System.out.println("Matching A Single Isolation Window ... WINDOW: [" + windowIndex + "] ...");

    ArrayList<Pair> isolationWindows = parseIsolationWindows();

    // Compute the MZ bounds of the Isolation Windows.
    double mzLower = (double) isolationWindows.get(windowIndex - 1).getL();
    double mzUpper = (double) isolationWindows.get(windowIndex - 1).getR();

    // TODO: This is Xiangyuan's Implementation, re-institute precursor searching
    String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
    String featureDetectParamFile = dataDirectory + "/featuredetect.params";
    String mzxmlFile = dataDirectory + "/r01_dia_data.mzXML";

    OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
    FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS2);
    featureDetect.init_MS2();


    /* Then */
    IsolationWindow isolationWindow = new IsolationWindow(mzLower, mzUpper); // test case for ArrayList<IsolationWindow>
    //featureDetect.detectMS2Features_method1("",isolationWindow.mzLow, isolationWindow.mzHigh, 20, 5e-6,5);
    int trail_rt_tolerance = 60; // 60s extension (theoretically should be infinite extension with this)
    int trail_rt_cutoff = 30; //30s hard cutoff
    double ms2_ppm_tolerance = 1e-5; // 10 ppm
    int min_peaknum_strict = 2; // imin 2 peaks
    double intensityNextPeakPercentageTol = 0.8;
    //featureDetect.detectMS2Features_method1(mzXMLFile,mzLower, mzUpper, trail_rt_tolerance, ms2_ppm_tolerance,5);
    boolean considerLessConfidentTrails = true;
    featureDetect.detectMS2Features_method2(mzxmlFile,mzLower,mzUpper,trail_rt_tolerance,ms2_ppm_tolerance,min_peaknum_strict,trail_rt_cutoff, intensityNextPeakPercentageTol, considerLessConfidentTrails);

//    featureDetect.detectMS2Features_method1(mzxmlFile, isolationWindow.mzLow, isolationWindow.mzHigh, 5, 5e-6, 4);
    isolationWindow.xics = featureDetect.getTrailsStrictTol();
    ArrayList<XIC> xicsLooseTol = featureDetect.getTrailsLooseTol();
    // Concatenate the Confident XICs to the non-confident ones
    isolationWindow.xics.addAll(xicsLooseTol);
    isolationWindow.removeDuplicateTrails();
    isolationWindow.xics.sort(Comparator.comparingDouble(XIC::getMZAtMaxIntensity));
    isolationWindow.populateMax();
    isolationWindow.populatertToTrailSet(0.16);
    isolationWindow.populatertToPrecursorSet(0.16, precursors);

    // will output to csv_data/ directory
    isolationWindow.writeIsolationWindowXICsToTSV(isolationWindow, this.WORKING_DIRECTORY, windowIndex);

    ArrayList<Peptide> candidatePeptides = this.parseFastaFile(this.MINIMUM_PEPTIDE_LENGTH);

    System.out.println("Digested Protein database to Peptides. Filtering peptides now ...");
    List<Peptide> peptideSubset =
        candidatePeptides.stream()
            .filter(
                x -> {
                  // Filter out ALL peptides that fall between the m/z of window bound
                  // Assuming we are just searching precursor charge 2 currently
                  double tmpMZ = Utils.MassToMz(x.mass, PRECURSOR_CHARGE);
                  return (mzLower < tmpMZ && mzUpper > tmpMZ);
                })
            .collect(Collectors.toList());

    // Export the Peptides
    RunFilter runFilter = new RunFilter();
    runFilter.writePeptidesToTSV(
        (ArrayList<Peptide>) peptideSubset,
        "/home/jia/Documents/Code/waterlooms/dia_data_reading/data/debug_output/");

    System.out.println("Completed Filtering, scoring now ...");

    int psetSize = peptideSubset.size();
    ArrayList<DbMatch> allMatches = new ArrayList<>(psetSize);

    for (int index = 0; index < peptideSubset.size(); index++) {
      if (index % 1000 == 0) {
        System.out.println("Matching Peptide: " + index + "/" + psetSize);
      }

      ArrayList<MatchedInterval> matchedIntervalsForPeptide =
          DbMatch.peptideMatch(peptideSubset.get(index), isolationWindow, isolationWindow, precursors); // TODO THIS IS A BUG SEE ABOVE

      // Compute the max score interval
        assert matchedIntervalsForPeptide != null;

        // TODO: uncomment this to write the debug_PEPTIDE_VALUE.tsv file in debug_output/
        //runFilter.writeMatchIntervalMatchedFeaturesTSV(matchedIntervalsForPeptide, this.WORKING_DIRECTORY);

      MatchedInterval maxScoreInterval = Collections.max(matchedIntervalsForPeptide, Comparator.comparing(x -> x.score));
        allMatches.add(
                new DbMatch(
                        peptideSubset.get(index).id,
                        peptideSubset.get(index).composition,
                        maxScoreInterval,peptideSubset.get(index))
        );
    }


    RetentionTimeMetadataTupleMap retentionTimeMetadataTupleMap = new RetentionTimeMetadataTupleMap();
    isolationWindow.trailToMatchedPeptidesSet.forEach(retentionTimeMetadataTupleMap::updateSetGivenTrail);
    retentionTimeMetadataTupleMap.normalizeRtTupleMapElements(isolationWindow.max);
    retentionTimeMetadataTupleMap.writeMapToTSV(windowIndex, this.WORKING_DIRECTORY);



    System.out.println("Completed Matching, sorting match subset now ...");

    Collections.sort(allMatches, Collections.reverseOrder());
    System.out.println("Sorting match completed, writing match subset now ...");
    System.out.println("Exiting before writing debug files.");
    System.exit(0);
    // Write allMatches to a file here for evaluation
    runFilter.writeUnscoredPeptidesToTSV(allMatches, this.WORKING_DIRECTORY);
    runFilter.writeExperimentResultTSV(allMatches,this.WORKING_DIRECTORY);
    runFilter.writeXICResultTSV(allMatches, this.WORKING_DIRECTORY);

    System.out.println("Done!");
  }

  private ArrayList<XIC> binarySearchPrecursorList(
      ArrayList<XIC> precursors, double mz, double rt) {

    final double PPM_TOLERANCE = 30;
    XIC upperMZXIC = new XIC(mz * (1 + PPM_TOLERANCE * Math.pow(10, -6)));
    XIC lowerMZXIC = new XIC(mz * (1 - PPM_TOLERANCE * Math.pow(10, -6)));

    // Heuristically cut the search space down just based on ppm
    int upperIndex =
        Collections.binarySearch(
            precursors, upperMZXIC, XIC::compareMzXIC);
    int lowerIndex =
        Collections.binarySearch(
            precursors, lowerMZXIC, XIC::compareMzXIC);
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

  private int parseInsertionPoint(int index) {
    if (index < 0) {
      return (-1 * (index) - 1);
    } else {
      return index;
    }
  }


  /**
   * Given an Absolute Path to a FeatureDetect Params File, return the Parameters in a Property
   * @param featureDetectParamFile String with absolute value of path
   * @return Properties of parameters file
   */
  private Properties parseFeatureDetectParamsFile(String featureDetectParamFile){
    Properties featureDetectParams = new Properties();
    try {
      FileReader fin = new FileReader(featureDetectParamFile);
      featureDetectParams.load(fin);
    } catch (IOException $e) {
      System.err.printf($e.getMessage());
      System.exit(1);
    }
    return featureDetectParams;
  }
}