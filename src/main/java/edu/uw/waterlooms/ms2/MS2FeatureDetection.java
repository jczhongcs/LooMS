package edu.uw.waterlooms.ms2;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.SVRScore;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.match.DbMatch;
import edu.uw.waterlooms.msutil.OpenMzxml;
import edu.uw.waterlooms.msutil.ScanEntry;
import edu.uw.waterlooms.msutil.MZXMLReader;
import java.io.*;
import java.util.*;

import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.util.Pair;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import umich.ms.fileio.filetypes.agilent.cef.jaxb.P;

public class MS2FeatureDetection {

  private static final int SECONDS_IN_MINUTES = 60;
  private static final int MS2_LEVEL = 2;

  public MS2FeatureDetection(
      int ms1Tolerance, String homePath, String mzXMLInFile, List<SVRScore> svrScores) {
    MS1TorPPM = ms1Tolerance;
    home_path = homePath;
    mzXMLFile = mzXMLInFile;
    precursorsFromSVR = svrScores;
  }

  // Overload constructor in the case that we only want MS2 XIC extraction
  public MS2FeatureDetection(int ms1Tolerance, String homePath, String mzXMLInFile) {
    MS1TorPPM = ms1Tolerance;
    home_path = homePath;
    mzXMLFile = mzXMLInFile;
  }

  static String home_path = "C:\\Users\\z46han\\Desktop\\Var_Temp\\Data\\";
  List<SVRScore> precursorsFromSVR;
  String mzXMLFile;
  // MS1 extraction
  private static int MS1TorPPM = 20;
  // global variables
  private static ArrayList<Double> mws = new ArrayList<>();
  private static ArrayList<double[]> dists = new ArrayList<>();

  private static double[] get_iso_dest(double mw, int num_iso) {
    int index = Collections.binarySearch(mws, mw);
    System.out.println("index: " + index);
    double[] mw_dist;
    // TODO: Possibly a bug but cap the index to be maximum of size 88 (thats how many elements are
    // in the averagine file)
    // TODO: Programmically calculate the distribution of the iso_dist at a given m/z instead of
    // using the averagine file

    if (index < 0) { // index = -(insertion point) - 1
      index = -(index + 1);
      index = (index) >= 89 ? 88 : index; // TODO: DEBUG THIS
      System.out.println("INDEX > 89: " + index);

      double[] iso_dist = Arrays.copyOfRange(dists.get(index), 0, num_iso);

      // int lowerIndexBound = (index - 1) < 0 ? 0 : (index - 1);

      int lowerIndexBound = Math.max(0, (index - 1));
      int upperIndexBound = Math.min(88, (index + 1));

      double[] iso_L = Arrays.copyOfRange(dists.get(lowerIndexBound), 0, num_iso);
      double[] iso_U = Arrays.copyOfRange(dists.get(upperIndexBound), 0, num_iso);
      double mw_L = mws.get(lowerIndexBound);
      double mw_U = mws.get(upperIndexBound);
      double[] x = {mw_L, mw, mw_U};
      mw_dist = new double[num_iso];
      for (int i = 0; i < num_iso; ++i) {
        double[] y = {iso_L[i], iso_dist[i], iso_U[i]};
        LinearInterpolator interp = new LinearInterpolator();
        PolynomialSplineFunction f = interp.interpolate(x, y);
        double dist_i = f.value(mw);
        mw_dist[i] = dist_i;
      }
    } else {
      index = (index) >= 89 ? 88 : index; // TODO: DEBUG THIS
      System.out.println("(negative insertion point) INDEX > 89: " + index);
      mw_dist = Arrays.copyOfRange(dists.get(index), 0, num_iso);
    }
    double dist_sum = Arrays.stream(mw_dist).sum();
    return Arrays.stream(mw_dist).map(x -> x / dist_sum).toArray();
  }

  private static int find_smallest(double[] a) {
    double smallest = MS1TorPPM;
    int index = -1;
    for (int i = 0; i < a.length; ++i) {
      if (a[i] <= smallest) {
        index = i;
        smallest = a[i];
      }
    }
    return index;
  }

  private static List<Double> extraction_ms1(
      List<IScan> spec_array, double[] iso_mass, double[] iso_dist) {
    int spec_size = spec_array.size();
    List<Double> extracted_ints = new ArrayList<>();
    for (int i = 0; i < spec_size; ++i) {
      extracted_ints.add(0.0);
    }
    for (int i = 0; i < iso_mass.length; ++i) {
      double mz = iso_mass[i];
      for (int j = 0; j < spec_size; ++j) {
        IScan spec = spec_array.get(j);
        ISpectrum spectrum = null;
        try {
          spectrum = spec.fetchSpectrum();
        } catch (FileParsingException e) {
          e.printStackTrace();
        }
        double[] mzs = spectrum.getMZs();
        double[] ints = spectrum.getIntensities();
        double[] mass_errors =
            Arrays.stream(mzs).map(s -> Math.abs((s - mz) / mz * Math.pow(10, 6))).toArray();
        int smallest_index = find_smallest(mass_errors);
        double extracted_int = smallest_index == -1 ? 0 : ints[smallest_index];
        extracted_int *= iso_dist[i];
        double cur_ints = extracted_ints.get(j);
        extracted_ints.set(j, cur_ints + extracted_int);
      }
    }
    return extracted_ints;
  }

  public void writeIsolationWindowPair(ArrayList<Pair<Double, Double>> isolationWindows, String pairOutfile) throws IOException {
    FileWriter fstream = new FileWriter(pairOutfile);
    BufferedWriter info = new BufferedWriter(fstream);
    int num_windows = isolationWindows.size();

    System.out.println("Writing IsolationWindow Ranges ...");

    for (int cur_iso = 0; cur_iso < num_windows; ++cur_iso) {
      info.write("START\n");
      info.write(
          Double.toString(isolationWindows.get(cur_iso).getKey())
              + " "
              + Double.toString(isolationWindows.get(cur_iso).getValue())
              + "\n");
      info.write("END\n");
    }
    info.close();
    fstream.close();
  }

  public IsolationWindow extractMS2XICForASingleIsolationWindow(int windowIndex, boolean OUTPUT_DEBUG) throws IOException {
    MZXMLReader mzml_reader = new MZXMLReader(mzXMLFile, MS2_LEVEL);
    ArrayList<Pair<Double, Double>> isolation_windows = mzml_reader.getIsolation_windows();
    ArrayList<ArrayList<ArrayList<ScanEntry>>> ms2_windowed = mzml_reader.getSpectrumPoints();

    ArrayList<ArrayList<ScanEntry>> isoWin = ms2_windowed.get(windowIndex);
    MS2Extraction ms2_extraction = new MS2Extraction(isoWin);
    Map<ScanEntry, XIC> MS2_xics = ms2_extraction.getMS2_xics();

    IsolationWindow isolationWindow = new IsolationWindow(
            isolation_windows.get(windowIndex).getKey(),
            isolation_windows.get(windowIndex).getValue()
    );
    isolationWindow.xics.addAll(MS2_xics.values());
    Collections.sort(isolationWindow.xics);
    isolationWindow.populateMax();


    if (OUTPUT_DEBUG){
      // Write values
    }
    return (isolationWindow);
  }


  /**
   * Extracts MS2 XICs from a given mzXML file.
   *
   * @throws IOException
   */
  public void extractMS2XIC(String workingDirectory, String outFile, int msLevel) throws Exception {
    MZXMLReader mzml_reader = new MZXMLReader(mzXMLFile, msLevel);
    ArrayList<Pair<Double, Double>> isolation_windows = mzml_reader.getIsolation_windows();
    ArrayList<ArrayList<ArrayList<ScanEntry>>> ms2_windowed = mzml_reader.getSpectrumPoints();
    System.out.println("Extracting MS2 XICs...");

    // For each IsolationWindow, write the MZ_LOW, MZ_HIGH for processing
    writeIsolationWindowPair(isolation_windows, outFile);

    // Export TRAILs detected for MS2
    int num_windows = isolation_windows.size();

    // Create the IsolationWindow Max Hashmap across ALL RT across ALL MZ
    HashMap<Double, Double> rtToMaxIntensity = new HashMap<>();

    for (int cur_iso = 0; cur_iso < num_windows; ++cur_iso) {
      System.out.println("Processing Isolation Window: [" + cur_iso + "/" + num_windows + "]");

      // Here is where you want to do the serialization code
      ArrayList<ArrayList<ScanEntry>> isoWin = ms2_windowed.get(cur_iso);

      // find local maximum ms2
      MS2Extraction ms2_extraction = new MS2Extraction(isoWin);
      // MS1 and MS2 XIC in this isolation window
      Map<ScanEntry, XIC> MS2_xics = ms2_extraction.getMS2_xics();

      IsolationWindow isolationWindow = new IsolationWindow(
          isolation_windows.get(cur_iso).getKey(),
          isolation_windows.get(cur_iso).getValue()
      );

      // foreach MS2_xics, append the value to the isolationWindow xics
      // TODO: This is a point of optimization, the XICs can be added as they're generated.
      isolationWindow.xics.addAll(MS2_xics.values());
      Collections.sort(isolationWindow.xics);

      isolationWindow.xics.forEach(
          xic -> {
            ArrayList rts = xic.getRetentionTimes();
            ArrayList its = xic.getIntensities();
            for (int rtIndex = 0; rtIndex < xic.getRetentionTimes().size(); rtIndex++){
              Double tmpRT = (double) rts.get(rtIndex);
              Double tmpIT = (double) its.get(rtIndex);

              // If key doesn't exist, add intensity : retentiontime map
              if (!rtToMaxIntensity.containsKey(tmpRT)){
                rtToMaxIntensity.put(tmpRT, tmpIT);
              } else {
                // If the key exists, replace it with a higher intensity one
                if (rtToMaxIntensity.get(tmpRT) < tmpIT)
                {
                  rtToMaxIntensity.replace(tmpRT, tmpIT);
                }
              }

            }

          }
      );

      // Call helper function to populate the max hashmap
      isolationWindow.populateMax();


      // Call Helper function to de-duplicate the trail(s) in a GREEDY way
      isolationWindow.removeDuplicateTrails();

      // Write Serialized IsolationWindow
      FileOutputStream isolationWindowOut = new FileOutputStream(
          workingDirectory + (cur_iso + 1) + "_window.ser"
      );
      ObjectOutputStream out = new ObjectOutputStream(isolationWindowOut);
      out.writeObject(isolationWindow);
      out.close();
      isolationWindowOut.close();

      // Foreach IsolationWindow, write the trail(s) to a CSV for R analysis & Matching
      writeIsolationWindowXICsToTSV(isolationWindow, workingDirectory, cur_iso);


    }

    // TODO: Depreciated serialization code
    /*
    // Write Serialized hashmap
    FileOutputStream retentionTimeToMaxIntensityOut = new FileOutputStream(
            workingDirectory + "retentionTimeToMaxIntensity.ser"
    );
    ObjectOutputStream outTwo = new ObjectOutputStream(retentionTimeToMaxIntensityOut);
    outTwo.writeObject(rtToMaxIntensity);
    outTwo.close();
    retentionTimeToMaxIntensityOut.close();
    System.out.println("Written the serialized map.");
    writeMaxHashmapToTSV(rtToMaxIntensity, workingDirectory);
    System.out.println("Written the tsv map.");
     */
  }

  private void writeMaxHashmapToTSV(
          HashMap<Double, Double> rtToMaxIntensityMap,
          String workingDirectory
  ) throws IOException {
    FileWriter fileWriter = new FileWriter(workingDirectory + "retentionTimeToMaxIntensity.tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

    bufferedWriter.write("retentionTime" + "\t" + "maxIntensity" + "\n");

    Iterator hmIterator = rtToMaxIntensityMap.entrySet().iterator();
    while (hmIterator.hasNext()){
      Map.Entry mapElement = (Map.Entry)hmIterator.next();
      bufferedWriter.write(mapElement.getKey() + "\t");
      bufferedWriter.write(mapElement.getValue() + "\n");
    }

    bufferedWriter.close();
    fileWriter.close();
  }

  private void writeIsolationWindowXICsToTSV(
      IsolationWindow isolationWindow, String workingDirectory, int windowIndex) throws IOException{
    FileWriter fileWriter = new FileWriter(workingDirectory + windowIndex + "_window.tsv");
    BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

    // Write the header
    // maxRT, maxMZ, maxINT, rts, mzs, ints
    bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\n");

    // For each XIC in XICs
    for (XIC xic : isolationWindow.xics) {
      // Write maxRT, maxMZ, maxINT
      bufferedWriter.write(Double.toString(xic.getRtAtMaxIntensity()) + "\t");
      bufferedWriter.write(Double.toString(xic.getMZAtMaxIntensity()) + "\t");
      bufferedWriter.write(Double.toString(xic.getMaxIntensity()) + "\t");

      // Write Retention Times
      StringBuilder stringBuilder = new StringBuilder();
      for (Double rts: xic.getRetentionTimes()){
        stringBuilder.append(Double.toString(rts * SECONDS_IN_MINUTES)).append(",");
      }
      bufferedWriter.write(stringBuilder.toString() + "\t");

      // Write MZs
      stringBuilder = new StringBuilder();
      for (Double mzs: xic.getMassChargeRatios()){
        stringBuilder.append(Double.toString(mzs)).append(",");
      }
      bufferedWriter.write(stringBuilder.toString() + "\t");

      // Write INTs
      stringBuilder = new StringBuilder();
      for (Double ints: xic.getIntensities()){
        stringBuilder.append(Double.toString(ints)).append(",");
      }
      bufferedWriter.write(stringBuilder.toString() + "\n");
    }

    bufferedWriter.close();
    fileWriter.close();
  }

  // Formerly featureDetect
  public void assemblePseudoSpectra(int msLevel) throws Exception {
    // read in mzml
    OpenMzxml of = new OpenMzxml(mzXMLFile);
    // Scans ms1_scans = new Scans(of, 1);
    MZXMLReader mzml_reader = new MZXMLReader(mzXMLFile, msLevel);
    ArrayList<Pair<Double, Double>> isolation_windows = mzml_reader.getIsolation_windows();
    int num_windows = isolation_windows.size();
    ArrayList<ArrayList<ArrayList<ScanEntry>>> ms2_windowed = mzml_reader.getSpectrumPoints();

    /*
    // TODO: Dynamically generate the averagine model as the pre-computed file is very limiting
    // TODO: The XIC is now done by Xiangyuan's Program, refactor to accept m/z, int, rt (1-4)
    // read in ms1 detected features
    CsvParserSettings settings = new CsvParserSettings();
    settings.detectFormatAutomatically();
    CsvParser parser = new CsvParser(settings);


    // MS1 extract XIC
    ArrayList<ArrayList<Map.Entry<Double, Precursor>>> ms1_xic_win = new ArrayList<>();
    for (int i = 0; i < num_windows; ++i) {
      ms1_xic_win.add(new ArrayList<>());
    }

    // Read in averagine (dynamically compute this in a later push)
    File file_avg = new File(home_path + "averagine_20160512.txt");
    BufferedReader br = new BufferedReader(new FileReader(file_avg));
    String line;
    while ((line = br.readLine()) != null) {
      String[] props = line.split("\t", 0);
      String mz = props[0];
      mws.add(Double.parseDouble(mz));
      String[] prop_dists = Arrays.copyOfRange(props, 1, props.length);
      double[] prop_dists_d = Arrays.stream(prop_dists).mapToDouble(Double::parseDouble).toArray();
      dists.add(prop_dists_d);
    }
    for (SVRScore precursor : precursorsFromSVR) {
      double p_mass_i = precursor.getMassChargeRatio(); // m/z
      int charge_i = (int) Math.round(precursor.getCharge()); // z
      double start_rt = precursor.getRetentionTimeStart(); // start RT
      double end_rt = precursor.getRetentionTimeEnd(); // end RT
      double peak_rt = precursor.getRetentionTime(); // rt
      double mw_i = charge_i * p_mass_i;
      double[] p_mass_iso = {p_mass_i};
      int num_iso = p_mass_iso.length;

      System.out.println("mw_i: " + mw_i + " num_iso: " + num_iso);
      // TODO: The ArrayOutOfBoundsError Is Occurring in get_iso_dest (Index: 89, Size: 89)
      // This will fail for very large precursors

      double[] ms1_dist = get_iso_dest(mw_i, num_iso);
      System.out.println("precursor: " + precursor.toString());
      // filter ms1 by rt
      Set<Map.Entry<Integer, IScan>> scans1 = ms1_scans.scanEntries;
      //System.out.println("scans1: " + scans1.toString());
      List<IScan> scans1_f =
          scans1.stream()
              .map(Map.Entry::getValue)
              .filter(s -> s.getRt() >= start_rt && s.getRt() <= end_rt)
              .collect(Collectors.toList());

      System.out.println("PRE-XIC: ");
      // get xic
      List<Double> rt =
          scans1.stream()
              .map(x -> x.getValue().getRt())
              .filter(s -> s >= start_rt && s <= end_rt)
              .collect(Collectors.toList());

      System.out.println("POST-XIC: ");

      List<Double> xic = extraction_ms1(scans1_f, p_mass_iso, ms1_dist);
      Precursor p = new Precursor(p_mass_i, charge_i, peak_rt, rt, xic);
      Map.Entry<Double, Precursor> xic_entry = new AbstractMap.SimpleEntry<>(peak_rt, p);
      // add to array list according to isolation window
      for (int j = 0; j < num_windows; ++j) {
        Pair<Double, Double> range = isolation_windows.get(j);
        // precursor in the j-th isolation window
        if (range.getKey() <= p_mass_i && p_mass_i <= range.getValue()) {
          ms1_xic_win.get(j).add(xic_entry);
          break;
        }
      }
    } */

    // Define the spectra.mgf output (not a mgf file, update the extension)
    FileWriter fstream = new FileWriter(home_path + "spectra.mgf");
    BufferedWriter info = new BufferedWriter(fstream);

    System.out.println("Extracting MS2 XICs...");
    // Export TRAILs detected for MS2
    for (int cur_iso = 0; cur_iso < num_windows; ++cur_iso) {
      System.out.println("Processing Isolation Window: [" + cur_iso + "/" + num_windows + "]");
      ArrayList<ArrayList<ScanEntry>> isoWin = ms2_windowed.get(cur_iso);
      // find local maximum ms2
      MS2Extraction ms2_extraction = new MS2Extraction(isoWin);
      // MS1 and MS2 XIC in this isolation window
      Map<ScanEntry, XIC> MS2_xics = ms2_extraction.getMS2_xics();

      info.write("START\n");
      // Write isolation window [start_mz, end_mz]
      info.write(
          Double.toString(isolation_windows.get(cur_iso).getKey())
              + " "
              + Double.toString(isolation_windows.get(cur_iso).getValue())
              + "\n");
      // For each ms2_extraction...
      Iterator xicIterator = MS2_xics.entrySet().iterator();
      while (xicIterator.hasNext()) {
        Map.Entry nextXIC = (Map.Entry) xicIterator.next();
        ScanEntry scanEntry = (ScanEntry) nextXIC.getKey();
        XIC xic = (XIC) nextXIC.getValue();

        // Write [local_max_mz, local_max_rt]
        info.write(
            Double.toString(scanEntry.mz)
                + " "
                + Double.toString(scanEntry.rt * SECONDS_IN_MINUTES)
                + "\n");
        // Write Trail Intensities
        StringBuilder builder = new StringBuilder();
        for (Double value : xic.getIntensities()) {
          builder.append(Double.toString(value)).append(" ");
        }
        info.write(builder.toString() + "\n");

        // Write Trail RTs
        builder = new StringBuilder();
        for (Double value : xic.getRetentionTimes()) {
          // Adjust the value by SECONDS_IN_MINUTES to get the RetentionTime in seconds
          builder.append(Double.toString(value * SECONDS_IN_MINUTES)).append(" ");
        }
        info.write(builder.toString() + "\n");
      }
      info.write("END\n");

      // TODO: Re-introduce this logic to generate pseudo-spectra
      // Currently PrecursorMatch is depreciated as we don't care about generating pseudo-spectra
      // new PrecursorMatch(MS2_xics, MS1_xics, info);
    }
    info.close();
  }
}
