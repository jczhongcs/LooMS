package edu.uw.waterlooms.ms1;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.Pair;
import edu.uw.waterlooms.entity.SignalIsotope;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.msutil.OpenMzxml;
import java.io.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

import edu.uw.waterlooms.peptideMatch.MSOneTrail;
import edu.uw.waterlooms.peptideMatch.MSOneTrailSet;
import edu.uw.waterlooms.peptideMatch.MSTwoTrail;
import edu.uw.waterlooms.peptideMatch.MSTwoTrailSet;
import edu.uw.waterlooms.service.ParameterService;
import org.apache.commons.lang3.time.StopWatch;
import org.json.JSONWriter;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import me.tongfei.progressbar.ProgressBar;

import javax.imageio.spi.IIOServiceProvider;

import static java.lang.Math.abs;
import static java.util.Arrays.binarySearch;

public class FeatureDetect {
  public enum DetectionType {
    MS1, MS2
  }
  private Set<Map.Entry<Integer, IScan>> scanEntries;
  private Set<Map.Entry<Integer, IScan>> scanEntries2;

  // MS1 params
  private double C12;
  private double C13;
  private double ppm;
  private double mzSearchRangePPM; // max tolerant m/z difference as to the same m/z
  private int MIN_CHECK; // number of consequent isotopic features to be identified
  private int MassStep;
  private double MinRelHeight;
  private int InvalidVal;
  private double PercentageOfMax;
  private double intensityThreshold;
  private int z_range;
  private int window_size;
  private int windowHi;
  private double IntensityBound;
  // MS1 features
  private double[][] intensityMap;
  private double[][] mzMap;
  private List<Integer>[] localpepMax;
  private int IDX;
  private double[] RT; // the array of real RT
  private double mzLow = Double.MAX_VALUE; // the smallest mz overall
  private double mzHi = Double.MIN_VALUE; // the largest mz overall
  private List<Double> scoreIntShape;
  private List<Double> scoreIsoDistr;
  private List<List<Pair<Double, Double>>> scorePepAll;
  private List<List<Pair<Integer, Integer>>>  scorePepAllPos;
  private List<Integer> scoreZVal;
  private List<Integer> scoreIsoNum;
  private List<List<Double>> intCluster;
  private List<Double> intSum;
  private List<Double> intensityPercentage;
  private List<List<XIC>> trailsAll;

  private int iSkipScanMS1;
//  private ArrayList<XIC> MS1Trails;

  // MS2 params
  // TODO: refactor parameters
  private int SEC_IN_MIN = 60;
  private int MZ_INDEX = 0;
  private int RT_INDEX = 1;
  private int INTENSITY_INDEX = 2;
  private int IFOCCUPIED_INDEX = 3;
  private double UNOCCUPIED = 0.0;
  private double OCCUPIED = 1.0;
  private double MZ_RANGE = 0.75; // for calculating intensityLocalAreaPercentage
  // MS2 features
  public HashMap<String, ArrayList<IScan>> scanByIsolationWindow;
//  private ArrayList<double[]> unassignedMS2Features;
  private List<double[]> unassignedMS2Features;
  private ArrayList<XIC> trailsStrictTol;
  private ArrayList<XIC> trailsLooseTol;
  /**
   * Read in all mzXML file information
   *
   * @param of new mzXML file opened by OpenMzxml
   */
  public FeatureDetect(OpenMzxml of, DetectionType detectionType) {
    switch (detectionType) {
      case MS1:
        setParameters();
        scanEntries = of.scanEntries;
        IDX = scanEntries.size();
        intensityMap = new double[IDX][];
        mzMap = new double[IDX][];
        RT = new double[IDX];

        // TODO: Added by jia to prevent init() from crashing
        scanEntries2 = of.scanEntries2;
        setiSkipScanMS1(of);
        break;
      case MS2:
        C12 = ParameterService.getC12Mass();
        C13 = ParameterService.getC13Mass();
        scanEntries = of.scanEntries;
        scanEntries2 = of.scanEntries2;
    }

    List<Integer> intKeysList = new ArrayList<>(of.num2scan2.keySet());
    System.out.println("intKeysList:"+JSONWriter.valueToString(intKeysList));
    System.out.println("iSkipScanMS1:"+iSkipScanMS1);

    System.out.println("scanEntries:"+JSONWriter.valueToString(of.num2scan2.get(intKeysList.get(0)).getPrecursor().getActivationInfo().getActivationEnergyLo()));
  }
  /** Helper function for constructor. Sets parameter values defined above. */
  private void setParameters() {
    // MS1 parameter
    C12 = ParameterService.getC12Mass();
    C13 = ParameterService.getC13Mass();
    ppm = ParameterService.getPPM();
    mzSearchRangePPM = ParameterService.getMzTolerancePpm();
    MIN_CHECK = ParameterService.getMinCheck();
    MassStep = ParameterService.getMassStep();
    MinRelHeight = ParameterService.getMinRelHeight();
    InvalidVal = ParameterService.getInvalidVal();
    PercentageOfMax = ParameterService.getPercentageOfMax();
    intensityThreshold = ParameterService.getIntensityThreshold();
    z_range = ParameterService.getzRange();
    window_size = ParameterService.getWindowSize();
    windowHi = ParameterService.getWindowHi();
    IntensityBound = ParameterService.getIntensityBound();
  }
  private void setiSkipScanMS1(OpenMzxml of ){
    if(of.num2scan!=null)
    {
      boolean bfirst = true;
      boolean bAllWindows = false;
      double dStratMZ=0.0;
      int iSkipScan= 0;
      int iWindowSize = 0;

      for (int i:of.num2scan.keySet()) {

        for (int scan2num : of.num2scan.get(i).getChildScans()) {
          double dMZlo = of.num2scan2.get(scan2num).getPrecursor().getMzRange().getLo();
          if(bfirst)
          {
            dStratMZ = dMZlo;
            bfirst = false;
          }
          else
          {
            if(dStratMZ==dMZlo)//已经循环了一次
            {
              bAllWindows =true;
            }
          }
          if(bAllWindows) {
            break;
          }
          iWindowSize++;
        }
        if (bAllWindows) {
          break;
        }
        iSkipScan++;
      }
      iSkipScanMS1 =  iSkipScan;
    }
  }



  public void init_MS2 () {
    scanByIsolationWindow = new HashMap<String, ArrayList<IScan>>();
    for (Map.Entry<Integer, IScan> scanEntry2 : scanEntries2) {
      IScan scan = scanEntry2.getValue();
      String windowRange = getWindowRange (scan.getPrecursor().getMzRangeStart(),  scan.getPrecursor().getMzRangeEnd());
      ArrayList<IScan> scanList;
      if (scanByIsolationWindow.containsKey(windowRange)) {
        scanList = scanByIsolationWindow.get(windowRange);
      } else {
        scanList = new ArrayList<>();
      }
      scanList.add(scan);
      scanByIsolationWindow.put(windowRange, scanList);
    }
  }
  private String getWindowRange (Double isoWindowLower, Double isoWindowUpper) {
    return "(" + isoWindowLower.intValue() + "," + isoWindowUpper.intValue() + ")";
  }

  /**
   * Method 2:
   *  Returns trailsStrictTol, a set of detected trails following strict constraints, and
   *     trailsLooseTol, a set of detected trails following loose constraints of the unused
   *     peaks from detecting trailsStrictTol.
   *
   * @param window
   * @param isolationWindowLower
   * @param isolationWindowUpper
   * @param rtNextPeakTolSec                adjacent peaks should have smaller rt than this value
   * @param mzTolerancePPMStrict            adjacent peaks should have smaller m/z than this value
   * @param MIN_PEAKNUM_Strict              total number of peaks in the trail should be larger than this value
   * @param rtMaxRangeSec                   length of the trail should not exceed this value; o/w trail is trimmied initialled from the local maximum
   * @param intensityNextPeakPercentageTol  allows for fluctuation of peak intensity in a trail
   */
  public void detectMS2Features_method2ForWindow(MSTwoTrailSet window, double isolationWindowLower, double isolationWindowUpper, double rtNextPeakTolSec, double mzTolerancePPMStrict,
                                                 int MIN_PEAKNUM_Strict, double rtMaxRangeSec, double intensityNextPeakPercentageTol, boolean considerLessConfidentTrails,
                                                 TreeMap<Double, List<MSOneTrail>> mapRTMSoneFeature){
    double mzTolerancePPMLoose = 10e-6; // 10ppm
    double  rtNextPeakTolMin = rtNextPeakTolSec/SEC_IN_MIN;
    double rtMaxRangeMin = rtMaxRangeSec/SEC_IN_MIN;

    ArrayList<IScan> scanList;
    String windowRange = getWindowRange(isolationWindowLower, isolationWindowUpper);
    scanList = scanByIsolationWindow.get(windowRange);
    // Check if the isolation window range is valid.
    if (!scanByIsolationWindow.containsKey(windowRange)) {
      System.err.printf("No isolation window range: " + windowRange);
      System.exit(1);
    }

    List<double[]> ms2InfoMap = Collections.synchronizedList(new ArrayList<>()); // Attempt to pre-alloc

    for (IScan scan: scanList) {
      double[] mzs = scan.getSpectrum().getMZs();
      double[] ints = scan.getSpectrum().getIntensities();
      double rt = scan.getRt();
      // release cache
      scan.setSpectrum(null, false);
      for (int p = 0; p < mzs.length; p++) {
        double[] ms2Info = new double[4];
        ms2Info[MZ_INDEX] = mzs[p];
        ms2Info[RT_INDEX] = rt;
        ms2Info[INTENSITY_INDEX] = ints[p];
        ms2Info[IFOCCUPIED_INDEX] = UNOCCUPIED;

        ms2InfoMap.add(ms2Info);
      }
    }
    System.gc();
    searchInfoList_method2ParallelForWindow(window,ms2InfoMap, rtNextPeakTolMin, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeMin, intensityNextPeakPercentageTol,mapRTMSoneFeature);
  }

  /**
   * Method 2:
   *  Returns trailsStrictTol, a set of detected trails following strict constraints, and
   *     trailsLooseTol, a set of detected trails following loose constraints of the unused
   *     peaks from detecting trailsStrictTol.
   *
   * @param filepath
   * @param isolationWindowLower
   * @param isolationWindowUpper
   * @param rtNextPeakTolSec                adjacent peaks should have smaller rt than this value
   * @param mzTolerancePPMStrict            adjacent peaks should have smaller m/z than this value
   * @param MIN_PEAKNUM_Strict              total number of peaks in the trail should be larger than this value
   * @param rtMaxRangeSec                   length of the trail should not exceed this value; o/w trail is trimmied initialled from the local maximum
   * @param intensityNextPeakPercentageTol  allows for fluctuation of peak intensity in a trail
   */
  public void detectMS2Features_method2(String filepath, double isolationWindowLower, double isolationWindowUpper, double rtNextPeakTolSec,  double mzTolerancePPMStrict, int MIN_PEAKNUM_Strict, double rtMaxRangeSec, double intensityNextPeakPercentageTol, boolean considerLessConfidentTrails){
    double mzTolerancePPMLoose = 10e-6; // 10ppm
    double  rtNextPeakTolMin = rtNextPeakTolSec/SEC_IN_MIN;
    double rtMaxRangeMin = rtMaxRangeSec/SEC_IN_MIN;

    ArrayList<IScan> scanList;
    String windowRange = getWindowRange(isolationWindowLower, isolationWindowUpper);
    scanList = scanByIsolationWindow.get(windowRange);
    // Check if the isolation window range is valid.
    if (!scanByIsolationWindow.containsKey(windowRange)) {
      System.err.printf("No isolation window range: " + windowRange);
      System.exit(1);
    }

    List<double[]> ms2InfoMap = Collections.synchronizedList(new ArrayList<>()); // Attempt to pre-alloc

    for (IScan scan: scanList) {
      double[] mzs = scan.getSpectrum().getMZs();
      double[] ints = scan.getSpectrum().getIntensities();
      double rt = scan.getRt();
      // release cache
      scan.setSpectrum(null, false);
      for (int p = 0; p < mzs.length; p++) {
        double[] ms2Info = new double[4];
        ms2Info[MZ_INDEX] = mzs[p];
        ms2Info[RT_INDEX] = rt;
        ms2Info[INTENSITY_INDEX] = ints[p];
        ms2Info[IFOCCUPIED_INDEX] = UNOCCUPIED;

        ms2InfoMap.add(ms2Info);
      }
    }
    System.gc();
    trailsStrictTol = searchInfoList_method2Parallel(ms2InfoMap, rtNextPeakTolMin, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeMin, intensityNextPeakPercentageTol);
//
////    unassignedMS2Features = new ArrayList<>();
//    unassignedMS2Features = Collections.synchronizedList(new ArrayList<>());
//    for (int i = 0; i < ms2InfoMap.size(); i++) {
//      if (ms2InfoMap.get(i)[IFOCCUPIED_INDEX] == UNOCCUPIED) {
//        unassignedMS2Features.add(ms2InfoMap.get(i));
//      }
//    }
//    if (considerLessConfidentTrails){
//      trailsLooseTol = searchInfoList_method2WithChargeParallel(unassignedMS2Features, rtNextPeakTolMin, mzTolerancePPMLoose, 1, rtMaxRangeMin, intensityNextPeakPercentageTol);
//    }
  }
  private ArrayList<XIC> searchInfoList_method2Parallel(List<double[]> infoList, double rtNextPeakTol, double mzTolerancePPM, int MIN_PEAKNUM_Strict,  double rtMaxRange, double intensityNextPeakPercentageTol) {
//    ProgressBar pb = new ProgressBar("Progress", infoList.size());
//    pb.start();

    ArrayList<XIC> trails = new ArrayList<>();
    sortInfoListParallel(infoList, MZ_INDEX);
    sortInfoListParallel(infoList, RT_INDEX);
    ArrayList<Pair<Integer, Double>> resultIndex; // A list of <index, rt> wrt the XIC.

    int long_trail_num = 0;
    int onepeak_trail_num = 0;

    int j;
    for (int i = 0; i < infoList.size();i++) { // Start the search from one peak (m/z smallest, this should not matter).
//      pb.step();

      if (infoList.get(i)[IFOCCUPIED_INDEX] == OCCUPIED) {
        continue; // Skip to the next if current is already in another group.
      }
      double prev_mz =  infoList.get(i)[MZ_INDEX];
      double prev_rt = infoList.get(i)[RT_INDEX];
      double prev_intensity = infoList.get(i)[INTENSITY_INDEX];
      boolean ifDecrease = false;
      double max_intensity = prev_intensity;
      double max_peak_rt = prev_rt;
      resultIndex = new ArrayList<>();
      resultIndex.add(new Pair<>(i, prev_rt));
      j = i + 1;
      while (j < infoList.size() && infoList.get(j)[IFOCCUPIED_INDEX] == OCCUPIED) {
        j++;
      }
      if (j == infoList.size()) {
        break;
      }

      for (; j < infoList.size(); j++) { // Continue the search from the next peak.
        double mz2 = infoList.get(j)[MZ_INDEX];
        double rt2 = infoList.get(j)[RT_INDEX];
        double intensity2 = infoList.get(j)[INTENSITY_INDEX];
        if (rt2 > prev_rt + rtNextPeakTol || rt2 > max_peak_rt + rtMaxRange / 2) {
          break;
        }
        if (infoList.get(j)[IFOCCUPIED_INDEX] == UNOCCUPIED
                && mz2 >= prev_mz * (1 - mzTolerancePPM)
                && mz2 <= prev_mz * (1 + mzTolerancePPM)
                && rt2 <= prev_rt + rtNextPeakTol
                && rt2 > prev_rt){
          if (intensity2 > prev_intensity && ifDecrease) {
            break;
          } else {
            ifDecrease = intensity2 < prev_intensity * intensityNextPeakPercentageTol;
          }
          prev_intensity = intensity2;
          resultIndex.add(new Pair<>(j, rt2));
          prev_mz = mz2;
          prev_rt = rt2;
          if (intensity2 > max_intensity) {
            max_intensity = intensity2;
            max_peak_rt = rt2;
          }
        }
      }
      // Make sure the left side is within rtMaxRange
      int k = 0;
      if (resultIndex.get(resultIndex.size()-1).getR() - resultIndex.get(0).getR() > rtMaxRange) {
        for (; k < resultIndex.size(); k++) {
          if (max_peak_rt - resultIndex.get(k).getR() < rtMaxRange / 2 ) {
            break;
          }
        }
      }
      // Must contain at least MIN_PEAKNUM_Strict peaks
      if (resultIndex.size() - k < MIN_PEAKNUM_Strict) {
        continue;
      }
      // Compute the number of long trails
      if (resultIndex.size() - k > 20) {
        long_trail_num++;
      }
      if (resultIndex.size() == 1) {
        onepeak_trail_num++;
      }

      // Compute the total intensity of the surrounding area of the trail: += MZ_RANGE
      double intensityLocalArea = 0;
      double benchmarkMZ = infoList.get(i)[MZ_INDEX];
      double rightRT = resultIndex.get(resultIndex.size()-1).getR();
      for (int p = i; p < infoList.size() && infoList.get(p)[RT_INDEX] <= rightRT; p++) {
        double curMZ = infoList.get(p)[MZ_INDEX];
        double curIntensity = infoList.get(p)[INTENSITY_INDEX];
        if (curMZ <= benchmarkMZ + MZ_RANGE && curMZ >= benchmarkMZ - MZ_RANGE) {
          intensityLocalArea += curIntensity;
        }
      }
      // Add XIC trail
      ArrayList<Double> mzs = new ArrayList<>();
      ArrayList<Double> rts = new ArrayList<>();
      ArrayList<Double> ints = new ArrayList<>();
      double totalIntensityOfTrail = 0;
      for (int x = k; x < resultIndex.size(); x++) { // Add to trails
        Pair<Integer, Double> item = resultIndex.get(x);
        int id = item.getL();
        mzs.add(infoList.get(id)[MZ_INDEX]);
        rts.add(infoList.get(id)[RT_INDEX]);
        ints.add(infoList.get(id)[INTENSITY_INDEX]);
        infoList.get(id)[IFOCCUPIED_INDEX] = OCCUPIED;

        totalIntensityOfTrail += infoList.get(id)[INTENSITY_INDEX];
      }

      //直接加载ms2TRAIL
      
      XIC trail = new XIC(ints, mzs, rts, totalIntensityOfTrail/intensityLocalArea);
      trail.quantifyPeaks();
      trails.add(trail);
    }
//    pb.stop();

    int count = 0;
    for (int i = 0; i < infoList.size();i++) {
      if (infoList.get(i)[IFOCCUPIED_INDEX] == OCCUPIED) {
        count++;
      }
    }
    double percentage = (double)count/infoList.size();
    System.out.println("Number of One peak trails = 1 peak: " + onepeak_trail_num);

    System.out.println("Number of long trails > 20 peaks: " + long_trail_num);
    System.out.println("Number of peaks covered: " + count + "/" + infoList.size());
    System.out.println("Percentage: " + percentage);
    System.out.println("Number of trails detected: " + trails.size());
    return trails;
  }
  private ArrayList<XIC> searchInfoList_method2ParallelForWindow(MSTwoTrailSet window, List<double[]> infoList, double rtNextPeakTol, double mzTolerancePPM,
                                                                 int MIN_PEAKNUM_Strict,  double rtMaxRange, double intensityNextPeakPercentageTol,
                                                                 TreeMap<Double, List<MSOneTrail>> mapRTMSoneFeature) {
//    ProgressBar pb = new ProgressBar("Progress", infoList.size());
//    pb.start();

    ArrayList<XIC> trails = new ArrayList<>();
    sortInfoListParallel(infoList, MZ_INDEX);
    sortInfoListParallel(infoList, RT_INDEX);
    ArrayList<Pair<Integer, Double>> resultIndex; // A list of <index, rt> wrt the XIC.

    int long_trail_num = 0;
    int onepeak_trail_num = 0;

    int j;
    for (int i = 0; i < infoList.size();i++) { // Start the search from one peak (m/z smallest, this should not matter).
//      pb.step();

      if (infoList.get(i)[IFOCCUPIED_INDEX] == OCCUPIED) {
        continue; // Skip to the next if current is already in another group.
      }
      double prev_mz =  infoList.get(i)[MZ_INDEX];
      double prev_rt = infoList.get(i)[RT_INDEX];
      double prev_intensity = infoList.get(i)[INTENSITY_INDEX];
      boolean ifDecrease = false;
      double max_intensity = prev_intensity;
      double max_peak_rt = prev_rt;
      resultIndex = new ArrayList<>();
      resultIndex.add(new Pair<>(i, prev_rt));
      j = i + 1;
      while (j < infoList.size() && infoList.get(j)[IFOCCUPIED_INDEX] == OCCUPIED) {
        j++;
      }
      if (j == infoList.size()) {
        break;
      }

      for (; j < infoList.size(); j++) { // Continue the search from the next peak.
        double mz2 = infoList.get(j)[MZ_INDEX];
        double rt2 = infoList.get(j)[RT_INDEX];
        double intensity2 = infoList.get(j)[INTENSITY_INDEX];
        if (rt2 > prev_rt + rtNextPeakTol || rt2 > max_peak_rt + rtMaxRange / 2) {
          break;
        }
        if (infoList.get(j)[IFOCCUPIED_INDEX] == UNOCCUPIED
                && mz2 >= prev_mz * (1 - mzTolerancePPM)
                && mz2 <= prev_mz * (1 + mzTolerancePPM)
                && rt2 <= prev_rt + rtNextPeakTol
                && rt2 > prev_rt){
          if (intensity2 > prev_intensity && ifDecrease) {
            break;
          } else {
            ifDecrease = intensity2 < prev_intensity * intensityNextPeakPercentageTol;
          }
          prev_intensity = intensity2;
          resultIndex.add(new Pair<>(j, rt2));
          prev_mz = mz2;
          prev_rt = rt2;
          if (intensity2 > max_intensity) {
            max_intensity = intensity2;
            max_peak_rt = rt2;
          }
        }
      }
      // Make sure the left side is within rtMaxRange
      int k = 0;
      if (resultIndex.get(resultIndex.size()-1).getR() - resultIndex.get(0).getR() > rtMaxRange) {
        for (; k < resultIndex.size(); k++) {
          if (max_peak_rt - resultIndex.get(k).getR() < rtMaxRange / 2 ) {
            break;
          }
        }
      }
      // Must contain at least MIN_PEAKNUM_Strict peaks
      if (resultIndex.size() - k < MIN_PEAKNUM_Strict) {
        continue;
      }
      // Compute the number of long trails
      if (resultIndex.size() - k > 20) {
        long_trail_num++;
      }
      if (resultIndex.size() == 1) {
        onepeak_trail_num++;
      }

      // Compute the total intensity of the surrounding area of the trail: += MZ_RANGE
      double intensityLocalArea = 0;
      double benchmarkMZ = infoList.get(i)[MZ_INDEX];
      double rightRT = resultIndex.get(resultIndex.size()-1).getR();
      for (int p = i; p < infoList.size() && infoList.get(p)[RT_INDEX] <= rightRT; p++) {
        double curMZ = infoList.get(p)[MZ_INDEX];
        double curIntensity = infoList.get(p)[INTENSITY_INDEX];
        if (curMZ <= benchmarkMZ + MZ_RANGE && curMZ >= benchmarkMZ - MZ_RANGE) {
          intensityLocalArea += curIntensity;
        }
      }
      // Add XIC trail
      double dints = 0.0;
      int ipos = 0;
      double[] mzs = new double[resultIndex.size()-k];
      double[] rts = new double[resultIndex.size()-k];
      double[] ints = new double[resultIndex.size()-k];
      double totalIntensityOfTrail = 0.0;
      double totalIntensityAreaOfTrail = 0.0;
      for (int x = k; x < resultIndex.size(); x++) { // Add to trails
        Pair<Integer, Double> item = resultIndex.get(x);
        int id = item.getL();
        mzs[x-k] = infoList.get(id)[MZ_INDEX];
        rts[x-k] = infoList.get(id)[RT_INDEX];
        ints[x-k] = infoList.get(id)[INTENSITY_INDEX];
        if(ints[x-k]>dints)
        {
          dints = ints[x-k];
          ipos = x-k;
        }
        if(x-k > 0) {
          totalIntensityAreaOfTrail += (ints[x-k] + ints[x-k-1]) * (rts[x-k] - rts[x-k -1]) / 2;
        }
        infoList.get(id)[IFOCCUPIED_INDEX] = OCCUPIED;

        totalIntensityOfTrail += infoList.get(id)[INTENSITY_INDEX];
      }

      //直接加载ms2TRAIL

//      MSTwoTrail trail = new MSTwoTrail(ints, mzs, rts, totalIntensityOfTrail/intensityLocalArea);
      window.readFromeTrail(mzs[ipos],rts[ipos],mzs,rts,ints,totalIntensityOfTrail,totalIntensityAreaOfTrail,mapRTMSoneFeature);
//      trail.quantifyPeaks();
//      trails.add(trail);
    }
    window.ms2Process();
//    pb.stop();

    int count = 0;
    for (int i = 0; i < infoList.size();i++) {
      if (infoList.get(i)[IFOCCUPIED_INDEX] == OCCUPIED) {
        count++;
      }
    }
    double percentage = (double)count/infoList.size();
    System.out.println("Number of One peak trails = 1 peak: " + onepeak_trail_num);

    System.out.println("Number of long trails > 20 peaks: " + long_trail_num);
    System.out.println("Number of peaks covered: " + count + "/" + infoList.size());
    System.out.println("Percentage: " + percentage);
    System.out.println("Number of trails detected: " + window.arrMSTwoTrail.size());
    return trails;
  }
  private ArrayList<XIC> searchInfoList_method2WithChargeParallel(List<double[]> infoList, double rtNextPeakTol, double mzTolerancePPM, int MIN_PEAKNUM_Strict,  double rtMaxRange, double intensityNextPeakPercentageTol) {
//    ProgressBar pb = new ProgressBar("Progress", infoList.size());
//    pb.start();


    ArrayList<XIC> trails = new ArrayList<>();
    sortInfoListParallel(infoList, MZ_INDEX);
    sortInfoListParallel(infoList, RT_INDEX);

      ArrayList<Pair<Integer, Double>> resultIndex; // A list of <index, rt> wrt the XIC.

    int long_trail_num = 0;
    int onepeak_trail_num = 0;

    int j;
    for (int i = 0; i < infoList.size();i++) { // Start the search from one peak (m/z smallest, this should not matter).

      if (infoList.get(i)[IFOCCUPIED_INDEX] == OCCUPIED) {
        continue; // Skip to the next if current is already in another group.
      }
      double prev_mz =  infoList.get(i)[MZ_INDEX];
      double dminMZ = (prev_mz * (1 - mzTolerancePPM) -  (C13 - C12));
      double dmaxMZ = (prev_mz * (1 + mzTolerancePPM) -  (C13 - C12));
      double prev_rt = infoList.get(i)[RT_INDEX];
      double prev_intensity = infoList.get(i)[INTENSITY_INDEX];
      boolean ifDecrease = false;
      double max_intensity = prev_intensity;
      double max_peak_rt = prev_rt;
      resultIndex = new ArrayList<>();
      resultIndex.add(new Pair<>(i, prev_rt));

      j = i + 1;
      while (j < infoList.size() && infoList.get(j)[IFOCCUPIED_INDEX] == OCCUPIED && infoList.get(j)[RT_INDEX] == prev_rt ) {
        j++;
      }
      if (j == infoList.size()) {
        break;
      }
      if (infoList.get(j)[RT_INDEX] > prev_rt) {
        continue;
      }

      for (; j < infoList.size() ; j++) { // Continue the search from the next peak.
        double mz2 = infoList.get(j)[MZ_INDEX];
        double rt2 = infoList.get(j)[RT_INDEX];
        double intensity2 = infoList.get(j)[INTENSITY_INDEX];
        if (rt2 > prev_rt ) {
          break;
        }
        if (infoList.get(j)[IFOCCUPIED_INDEX] == UNOCCUPIED
                && mz2 >= dminMZ
                && mz2 <= dmaxMZ
                && rt2 == prev_rt ){
//          if (intensity2 < prev_intensity) {
          if (
                  (prev_mz<1500 &&   intensity2 < prev_intensity)
                          ||
                          (prev_mz>=1500 && intensity2/2.0  <  prev_intensity)
          ) {
            break;
          }
//          System.out.println(mz2+"\t"+prev_mz+"\t"+dminMZ+"\t"+dmaxMZ+"\t"+rt2+"\t"+prev_rt);
          prev_intensity = intensity2;
          resultIndex.add(new Pair<>(j, rt2));
          prev_mz = mz2;
          prev_rt = rt2;
          dminMZ = (prev_mz * (1 - mzTolerancePPM) -  (C13 - C12));
          dmaxMZ = (prev_mz * (1 + mzTolerancePPM) -  (C13 - C12));
          if (intensity2 > max_intensity) {
            max_intensity = intensity2;
            max_peak_rt = rt2;
          }
        }
      }
      // Make sure the left side is within rtMaxRange
      int k = 0;
      if (resultIndex.size() - k > 20) {
        long_trail_num++;
      }
      if (resultIndex.size() == 1) {
        continue;
//        onepeak_trail_num++;
      }

      // Compute the total intensity of the surrounding area of the trail: += MZ_RANGE
      double intensityLocalArea = 0;
      // Add XIC trail
      ArrayList<Double> mzs = new ArrayList<>();
      ArrayList<Double> rts = new ArrayList<>();
      ArrayList<Double> ints = new ArrayList<>();
      double totalIntensityOfTrail = 0;
      for (int x = k; x < resultIndex.size(); x++) { // Add to trails
        Pair<Integer, Double> item = resultIndex.get(x);
        int id = item.getL();
        infoList.get(id)[IFOCCUPIED_INDEX] = OCCUPIED;

      }
      Pair<Integer, Double> item = resultIndex.get(resultIndex.size()-1);
      int firstid = item.getL();

      mzs.add(infoList.get(firstid)[MZ_INDEX]);
      rts.add(infoList.get(firstid)[RT_INDEX]);
      ints.add(infoList.get(firstid)[INTENSITY_INDEX]);

      totalIntensityOfTrail = infoList.get(firstid)[INTENSITY_INDEX];
      XIC trail = new XIC(ints, mzs, rts, totalIntensityOfTrail);
      trail.quantifyPeaks();
      trails.add(trail);
    }
    int count = trails.size();
    double percentage = (double)count/infoList.size();

    System.out.println("Number of long trails > 20 peaks: " + long_trail_num);
    System.out.println("Number of peaks covered: " + count + "/" + infoList.size());
    System.out.println("Percentage: " + percentage);
    System.out.println("Number of trails detected: " + trails.size());
    return trails;
  }

  private ArrayList<XIC> searchInfoList_method2(ArrayList<double[]> infoList, double rtNextPeakTol, double mzTolerancePPM, int MIN_PEAKNUM_Strict,  double rtMaxRange, double intensityNextPeakPercentageTol) {
//    ProgressBar pb = new ProgressBar("Progress", infoList.size());
//    pb.start();

    ArrayList<XIC> trails = new ArrayList<>();
    sortInfoList(infoList, MZ_INDEX);
    sortInfoList(infoList, RT_INDEX);
    ArrayList<Pair<Integer, Double>> resultIndex; // A list of <index, rt> wrt the XIC.

    int long_trail_num = 0;
    int onepeak_trail_num = 0;

    int j;
    for (int i = 0; i < infoList.size();i++) { // Start the search from one peak (m/z smallest, this should not matter).
//      pb.step();

      if (infoList.get(i)[IFOCCUPIED_INDEX] == OCCUPIED) {
        continue; // Skip to the next if current is already in another group.
      }
      double prev_mz =  infoList.get(i)[MZ_INDEX];
      double prev_rt = infoList.get(i)[RT_INDEX];
      double prev_intensity = infoList.get(i)[INTENSITY_INDEX];
      boolean ifDecrease = false;
      double max_intensity = prev_intensity;
      double max_peak_rt = prev_rt;
      resultIndex = new ArrayList<>();
      resultIndex.add(new Pair<>(i, prev_rt));
      j = i + 1;
      while (j < infoList.size() && infoList.get(j)[IFOCCUPIED_INDEX] == OCCUPIED) {
        j++;
      }
      if (j == infoList.size()) {
        break;
      }

      for (; j < infoList.size(); j++) { // Continue the search from the next peak.
        double mz2 = infoList.get(j)[MZ_INDEX];
        double rt2 = infoList.get(j)[RT_INDEX];
        double intensity2 = infoList.get(j)[INTENSITY_INDEX];
        if (rt2 > prev_rt + rtNextPeakTol || rt2 > max_peak_rt + rtMaxRange / 2) {
          break;
        }
        if (infoList.get(j)[IFOCCUPIED_INDEX] == UNOCCUPIED
                && mz2 >= prev_mz * (1 - mzTolerancePPM)
                && mz2 <= prev_mz * (1 + mzTolerancePPM)
                && rt2 <= prev_rt + rtNextPeakTol
                && rt2 > prev_rt){
          if (intensity2 > prev_intensity && ifDecrease) {
            break;
          } else {
            ifDecrease = intensity2 < prev_intensity * intensityNextPeakPercentageTol;
          }
          prev_intensity = intensity2;
          resultIndex.add(new Pair<>(j, rt2));
          prev_mz = mz2;
          prev_rt = rt2;
          if (intensity2 > max_intensity) {
            max_intensity = intensity2;
            max_peak_rt = rt2;
          }
        }
      }
      // Make sure the left side is within rtMaxRange
      int k = 0;
      if (resultIndex.get(resultIndex.size()-1).getR() - resultIndex.get(0).getR() > rtMaxRange) {
        for (; k < resultIndex.size(); k++) {
          if (max_peak_rt - resultIndex.get(k).getR() < rtMaxRange / 2 ) {
            break;
          }
        }
      }
      // Must contain at least MIN_PEAKNUM_Strict peaks
      if (resultIndex.size() - k < MIN_PEAKNUM_Strict) {
        continue;
      }
      // Compute the number of long trails
      if (resultIndex.size() - k > 20) {
        long_trail_num++;
      }
      if (resultIndex.size() == 1) {
        onepeak_trail_num++;
      }

      // Compute the total intensity of the surrounding area of the trail: += MZ_RANGE
      double intensityLocalArea = 0;
      double benchmarkMZ = infoList.get(i)[MZ_INDEX];
      double rightRT = resultIndex.get(resultIndex.size()-1).getR();
      for (int p = i; p < infoList.size() && infoList.get(p)[RT_INDEX] <= rightRT; p++) {
        double curMZ = infoList.get(p)[MZ_INDEX];
        double curIntensity = infoList.get(p)[INTENSITY_INDEX];
        if (curMZ <= benchmarkMZ + MZ_RANGE && curMZ >= benchmarkMZ - MZ_RANGE) {
          intensityLocalArea += curIntensity;
        }
      }
      // Add XIC trail
      ArrayList<Double> mzs = new ArrayList<>();
      ArrayList<Double> rts = new ArrayList<>();
      ArrayList<Double> ints = new ArrayList<>();
      double totalIntensityOfTrail = 0;
      for (int x = k; x < resultIndex.size(); x++) { // Add to trails
        Pair<Integer, Double> item = resultIndex.get(x);
        int id = item.getL();
        mzs.add(infoList.get(id)[MZ_INDEX]);
        rts.add(infoList.get(id)[RT_INDEX]);
        ints.add(infoList.get(id)[INTENSITY_INDEX]);
        infoList.get(id)[IFOCCUPIED_INDEX] = OCCUPIED;

        totalIntensityOfTrail += infoList.get(id)[INTENSITY_INDEX];
      }
      XIC trail = new XIC(ints, mzs, rts, totalIntensityOfTrail/intensityLocalArea);
      trail.quantifyPeaks();
      trails.add(trail);
    }
//    pb.stop();

    int count = 0;
    for (int i = 0; i < infoList.size();i++) {
      if (infoList.get(i)[IFOCCUPIED_INDEX] == OCCUPIED) {
        count++;
      }
    }
    double percentage = (double)count/infoList.size();
    System.out.println("Number of One peak trails = 1 peak: " + onepeak_trail_num);

    System.out.println("Number of long trails > 20 peaks: " + long_trail_num);
    System.out.println("Number of peaks covered: " + count + "/" + infoList.size());
    System.out.println("Percentage: " + percentage);
    System.out.println("Number of trails detected: " + trails.size());
    return trails;
  }

  public ArrayList<XIC> getTrailsStrictTol() {
    return trailsStrictTol;
  }

  public ArrayList<XIC> getTrailsLooseTol() {
    return trailsLooseTol;
  }

  /**
   * init() not needed
   *
   * @param rtNextPeakTolSec
   * @param mzTolerancePPMStrict
   * @param MIN_PEAKNUM_Strict
   * @param rtMaxRangeSec
   * @param intensityNextPeakPercentageTol
   */
  public void detectMS1Trails(String mzxmlFile, double rtNextPeakTolSec,  double mzTolerancePPMStrict, int MIN_PEAKNUM_Strict, double rtMaxRangeSec, double intensityNextPeakPercentageTol){
    double  rtNextPeakTolMin = rtNextPeakTolSec/SEC_IN_MIN;
    double rtMaxRangeMin = rtMaxRangeSec/SEC_IN_MIN;

    int numMzs = 0;
    ArrayList<double[]> ms1InfoMap = new ArrayList<>();
    for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
      IScan scan = scanEntry.getValue();
      ISpectrum spectrum = null;
      try {
        spectrum = scan.fetchSpectrum();
      } catch (FileParsingException e) {
        e.printStackTrace();
      }
      double rt = scan.getRt();
      double[] mzs = spectrum.getMZs();
      double[] ints = spectrum.getIntensities();
      numMzs += mzs.length;

      for (int p = 0; p < mzs.length; p++) {
        double[] ms1Info = new double[4];
        ms1Info[MZ_INDEX] = mzs[p];
        ms1Info[RT_INDEX] = rt;
        ms1Info[INTENSITY_INDEX] = ints[p];
        ms1Info[IFOCCUPIED_INDEX] = UNOCCUPIED;

        ms1InfoMap.add(ms1Info);
      }
    }

    System.out.println("NumScans & NumMZs " + scanEntries.size() + " | " + numMzs);
//    MS1Trails = searchInfoList_method2(ms1InfoMap, rtNextPeakTolMin, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeMin, intensityNextPeakPercentageTol);
  }

//  public ArrayList<XIC> getMS1Trails() {
//    return this.MS1Trails;
//  }

  public void writeMS1TrailsData(String filename, ArrayList<XIC> trails) throws IOException {
    File file = new File(filename);
    try {
      // create FileWriter object with file as parameter
      FileWriter outputfile = new FileWriter(file);
      // PrintWriter
      PrintWriter printWriter = new PrintWriter(outputfile);
      String header = "mz\trt\tpeak_sum\tpeaks_area\tints\n";
      printWriter.print(header);

      for (int i = 0; i < trails.size(); i++) {
        StringBuilder data = new StringBuilder();
        XIC trail = trails.get(i);
        trail.quantifyPeaks();
        data.append(trail.getMZAtMaxIntensity()).append('\t');
        data.append(trail.getRtAtMaxIntensity()).append('\t');
        data.append(trail.getPeakSum()).append('\t');
        data.append(trail.getPeakArea()).append('\t');
        data.append('\n');
        printWriter.print(data);
      }
      // closing writer connection
      printWriter.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public void writeMS2TrailsData(String filename, ArrayList<XIC> trails) throws IOException {
        File file = new File(filename);
        try {
            // create FileWriter object with file as parameter
            FileWriter outputfile = new FileWriter(file);
            // PrintWriter
            PrintWriter printWriter = new PrintWriter(outputfile);
            String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";
            printWriter.print(header);

            for (int i = 0; i < trails.size(); i++) {
                StringBuilder data = new StringBuilder();
                ArrayList<Double> mzs =  trails.get(i).getMassChargeRatios();
                ArrayList<Double> rts = trails.get(i).getRetentionTimes();
                ArrayList<Double> ints = trails.get(i).getIntensities();

                data.append(trails.get(i).getMZAtMaxIntensity()).append('\t');
                data.append(trails.get(i).getRtAtMaxIntensity()).append('\t');
                String mzsStr = mzs.get(0).toString();
                String rtsStr = rts.get(0).toString();
                String intsStr = ints.get(0).toString();
                double peaks_sum = ints.get(0);;
                double peaks_area = 0;
                for (int j = 1; j < ints.size(); j++) {
                    peaks_sum += ints.get(j);
                    peaks_area += (ints.get(j) + ints.get(j-1)) * (rts.get(j) - rts.get(j-1)) / 2;
                    mzsStr += ',' + mzs.get(j).toString();
                    rtsStr += ',' + rts.get(j).toString();
                    intsStr += ',' + ints.get(j).toString();
                }
                data.append(mzsStr).append('\t');
                data.append(rtsStr).append('\t');
                data.append(intsStr).append('\t');
                data.append(peaks_sum).append('\t');
                data.append(peaks_area).append('\t');

                data.append('\n');
                printWriter.print(data);
            }
            // closing writer connection
          printWriter.flush();
          printWriter.close();
          outputfile.close();
//            printWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


  public void writeMS2TrailsData(String filename, ArrayList<XIC> trails,ArrayList<XIC> trailsLoose) throws IOException {
    File file = new File(filename);
    double drtspan = -1.0;
    try {
      // create FileWriter object with file as parameter
      FileWriter outputfile = new FileWriter(file);
      // PrintWriter
      PrintWriter printWriter = new PrintWriter(outputfile);
      String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";
      printWriter.print(header);

      for (int i = 0; i < trails.size(); i++) {
        StringBuilder data = new StringBuilder();
        ArrayList<Double> mzs =  trails.get(i).getMassChargeRatios();
        ArrayList<Double> rts = trails.get(i).getRetentionTimes();
        ArrayList<Double> ints = trails.get(i).getIntensities();

        data.append(trails.get(i).getMZAtMaxIntensity()).append('\t');
        data.append(trails.get(i).getRtAtMaxIntensity()).append('\t');
        String mzsStr = mzs.get(0).toString();
        String rtsStr = rts.get(0).toString();
        String intsStr = ints.get(0).toString();
        double peaks_sum = ints.get(0);;
        double peaks_area = 0;
        for (int j = 1; j < ints.size(); j++) {
          if(drtspan<0) {
            drtspan = rts.get(j) - rts.get(j - 1);
          }
          peaks_sum += ints.get(j);
          peaks_area += (ints.get(j) + ints.get(j-1)) * (rts.get(j) - rts.get(j-1)) / 2;
          mzsStr += ',' + mzs.get(j).toString();
          rtsStr += ',' + rts.get(j).toString();
          intsStr += ',' + ints.get(j).toString();
        }
        data.append(mzsStr).append('\t');
        data.append(rtsStr).append('\t');
        data.append(intsStr).append('\t');
        data.append(peaks_sum).append('\t');
        data.append(peaks_area).append('\t');

        data.append('\n');
        printWriter.print(data);
      }
      if(trailsLoose!=null)
      {
        for (int i = 0; i < trailsLoose.size(); i++) {
          StringBuilder data = new StringBuilder();
          ArrayList<Double> mzs =  trailsLoose.get(i).getMassChargeRatios();
          ArrayList<Double> rts = trailsLoose.get(i).getRetentionTimes();
          ArrayList<Double> ints = trailsLoose.get(i).getIntensities();

          data.append(trailsLoose.get(i).getMZAtMaxIntensity()).append('\t');
          data.append(trailsLoose.get(i).getRtAtMaxIntensity()).append('\t');
          String mzsStr = mzs.get(0).toString();
          String rtsStr = rts.get(0).toString();
          String intsStr = ints.get(0).toString();
          double peaks_sum = ints.get(0);;
          double peaks_area = 0;
          for (int j = 1; j < ints.size(); j++) {
            peaks_sum += ints.get(j);
            peaks_area += (ints.get(j) + ints.get(j-1)) * (rts.get(j) - rts.get(j-1)) / 2;
            mzsStr += ',' + mzs.get(j).toString();
            rtsStr += ',' + rts.get(j).toString();
            intsStr += ',' + ints.get(j).toString();
          }
          //just a peak
          if(ints.size()==1) {
            peaks_area += (peaks_sum * drtspan) / 2;
          }

          data.append(mzsStr).append('\t');
          data.append(rtsStr).append('\t');
          data.append(intsStr).append('\t');
          data.append(peaks_sum).append('\t');
          data.append(peaks_area).append('\t');

          data.append('\n');
          printWriter.print(data);
        }
      }

      // closing writer connection
      printWriter.flush();
      printWriter.close();
      outputfile.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
    /**
   * Call this function to detect features from LC-MS, write to a file, followed by machine learning
   * scoring methods
   */
  public void detectFeatures(String oldFilePath) {
      String filepath = oldFilePath.replaceFirst("[.][^.]+$", "");
    StopWatch ms1FeatureinitStopWatch = new StopWatch();
    ms1FeatureinitStopWatch.start();
    init_MS1();

    ms1FeatureinitStopWatch.stop();
    System.out.println("Init  mzML File Elapsed Time in Minutes: " + ms1FeatureinitStopWatch.getTime(TimeUnit.MINUTES) + " ...");

//      findLocalMax();
    StopWatch ms1FeaturefindMaxStopWatch = new StopWatch();
    ms1FeaturefindMaxStopWatch.start();
      findLocalMaxWithSkipScan();
    ms1FeaturefindMaxStopWatch.stop();
    System.out.println("Findmax  mzML File Elapsed Time in Minutes: " + ms1FeaturefindMaxStopWatch.getTime(TimeUnit.MINUTES) + " ...");

    System.out.println("findLocalMax completed");
//      searchIsotope(window_size, z_range);


      searchIsotopeWithSkipScan(window_size, z_range);


      System.out.println("searchIsotope completed");
      isoDistributionScore();
      System.out.println("isoDistributionScore completed");
//      detectTrailofPrecursor();
      detectTrailofPrecursorWithSkipScan();
      System.out.println("detectTrailofPrecursor completed");
      try {
        writeFeatures(
                filepath
                        + "_featureAllZ.tsv"); // since a feature is searched with all possible charge states
      } catch (IOException e) {
        System.out.println("Write to file failed.");
      }
    }

  /** Read in for each Spectrum, including intensities of recorded m/z values */
//  private void init_MS1() {
  public void init_MS1() {
    int count = 0;

    double rtHi = 0;
    for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
      IScan scan = scanEntry.getValue();
      ISpectrum spectrum = null;
      try {
        spectrum = scan.fetchSpectrum();
      } catch (FileParsingException e) {
        e.printStackTrace();
      }
      double[] mz = spectrum.getMZs();
      double[] intensities = spectrum.getIntensities();
      if (mz[mz.length - 1] > mzHi) {
        mzHi = mz[mz.length - 1];
      }
      if (mz[0] < mzLow) {
        mzLow = mz[0];
      }
      RT[count] = scan.getRt();
      intensityMap[count] = intensities;
      mzMap[count] = mz;
      count++;

      if (rtHi < scan.getRt()) {
        rtHi = scan.getRt();
      }
    }
    // Choose window size; must be in 3-15 scans
    int windowCalculated = (int) Math.ceil(rtHi/15); // scaling function
    if (windowCalculated > window_size) {
      window_size = windowCalculated;
    }
    if (window_size > windowHi) {
      window_size = windowHi;
    }
    for (Map.Entry<Integer, IScan> scanEntry2 : scanEntries2) {
      IScan scan = scanEntry2.getValue();
      Object scanPrecursorCharge = scan.getPrecursor().getCharge();
      if (scanPrecursorCharge != null && z_range < (int)scanPrecursorCharge) {
        z_range = scan.getPrecursor().getCharge();
      }else {
        z_range=8;
      }
    }
    System.out.println("window_size: " + window_size);
    System.out.println("z_range: " + z_range);
  }


  /** Initial detection on local maximal intensities on features */
  private void findLocalMaxWithSkipScan(){
    // Copy of intensityMap
    // Records if a m/z is reached. If reached, set to -1
    double[][] intMap = new double[IDX][];
    for (int i = 0; i < IDX; i++) {
      double[] copy = new double[intensityMap[i].length];
      for (int j = 0; j < intensityMap[i].length; j++) {
        copy[j] = intensityMap[i][j];
      }
      intMap[i] = copy;
    }

    localpepMax = new List[IDX];
    for(int i = 0; i < IDX; i++){
      localpepMax[i] = new ArrayList<>(); // initiate local max map
    }

    // Search max near one m/z horizontally
//    System.out.println("idx = " + IDX);
    for(int width = 0; width < IDX; width++){
      for(int height = 0; height < intMap[width].length; height++){
        if (intMap[width][height] == InvalidVal) {
          continue;
        }
        int count = 1;
        int col_for_max = width; // Column for max
        int pos_for_max = height; // Position for max
        int prev_col_for_max = InvalidVal; // Column for max
        int prev_pos_for_max = InvalidVal; // Position for max
        double pit_intval = InvalidVal;
        int last_col = width;
        int last_pos = height;
        double value_int_max = intMap[col_for_max][pos_for_max]; // The max value of intensity of one peptide
        double prev_value_int_max = InvalidVal;
        double prev_val = intMap[last_col][last_pos];
        intMap[col_for_max][pos_for_max] = InvalidVal; // Change the visited position to invalid

        // Search at the right column
        int rightPos = 0;
        boolean isDecrease = false;
//        for (int rightCol = width + 1; rightCol < IDX; rightCol++){
        for (int rightCol = width + iSkipScanMS1; rightCol < IDX; rightCol+=iSkipScanMS1){
          rightPos = searchRange(mzMap[rightCol], mzMap[last_col][last_pos] * (1-mzSearchRangePPM), mzMap[last_col][last_pos] * (1+mzSearchRangePPM));
          // multi-local-max
          //start
          if (rightPos == InvalidVal){
            if (count >= MIN_CHECK) {
              localpepMax[col_for_max].add(pos_for_max);
              if (prev_value_int_max!=InvalidVal && pit_intval < value_int_max * PercentageOfMax) {
                localpepMax[prev_col_for_max].add(prev_pos_for_max);
              }
            }
            count = 1;
            break;
          }
          else if (intMap[rightCol][rightPos] < prev_val) {
            if (value_int_max <= prev_value_int_max) {
              value_int_max = prev_value_int_max;
              col_for_max = prev_col_for_max;
              pos_for_max = prev_pos_for_max;
            }
            else if (prev_value_int_max!=InvalidVal && pit_intval < value_int_max * PercentageOfMax) { // the second max, value_int_max > prev_value_int_max
              localpepMax[prev_col_for_max].add(prev_pos_for_max);
            }
            prev_value_int_max = InvalidVal; // turned off when not necessary
            isDecrease = true;
          }
          else if (isDecrease) { // cur_val >= prev_val
            if (prev_val < value_int_max * PercentageOfMax) { // prev_val is lowest intensity point
              localpepMax[col_for_max].add(pos_for_max);
            } else {
              prev_col_for_max = col_for_max;
              prev_pos_for_max = pos_for_max;
              prev_value_int_max = value_int_max;
            }
            value_int_max = intMap[rightCol][rightPos];
            col_for_max = rightCol;
            pos_for_max = rightPos;
            pit_intval = prev_val;
            isDecrease = false;
          }
          if (intMap[rightCol][rightPos] > value_int_max) {
            value_int_max = intMap[rightCol][rightPos];
            col_for_max = rightCol;
            pos_for_max = rightPos;
          }
          // update
          last_col = rightCol;
          last_pos = rightPos;
          prev_val = intMap[rightCol][rightPos];
          count++;
          //end

          // turn the visited to invalid
          intMap[rightCol][rightPos] = InvalidVal;

        }
        // add last local max position;
        if (rightPos == IDX && count >= MIN_CHECK) {
          localpepMax[col_for_max].add(pos_for_max);
        }
      }
    }
    // sort
    for (int i = 0; i < IDX; i++) {
      Collections.sort(localpepMax[i]);
    }
  }
  /** Initial detection on local maximal intensities on features */
  private void findLocalMax(){
    // Copy of intensityMap
    // Records if a m/z is reached. If reached, set to -1
    double[][] intMap = new double[IDX][];
    for (int i = 0; i < IDX; i++) {
      double[] copy = new double[intensityMap[i].length];
      for (int j = 0; j < intensityMap[i].length; j++) {
        copy[j] = intensityMap[i][j];
      }
      intMap[i] = copy;
    }

    localpepMax = new List[IDX];
    for(int i = 0; i < IDX; i++){
      localpepMax[i] = new ArrayList<>(); // initiate local max map
    }

    // Search max near one m/z horizontally
//    System.out.println("idx = " + IDX);
    for(int width = 0; width < IDX; width++){
      for(int height = 0; height < intMap[width].length; height++){
        if (intMap[width][height] == InvalidVal) {
          continue;
        }
        int count = 1;
        int col_for_max = width; // Column for max
        int pos_for_max = height; // Position for max
        int prev_col_for_max = InvalidVal; // Column for max
        int prev_pos_for_max = InvalidVal; // Position for max
        double pit_intval = InvalidVal;
        int last_col = width;
        int last_pos = height;
        double value_int_max = intMap[col_for_max][pos_for_max]; // The max value of intensity of one peptide
        double prev_value_int_max = InvalidVal;
        double prev_val = intMap[last_col][last_pos];
        intMap[col_for_max][pos_for_max] = InvalidVal; // Change the visited position to invalid

        // Search at the right column
        int rightPos = 0;
        boolean isDecrease = false;
        for (int rightCol = width + 1; rightCol < IDX; rightCol++){
          rightPos = searchRange(mzMap[rightCol], mzMap[last_col][last_pos] * (1-mzSearchRangePPM), mzMap[last_col][last_pos] * (1+mzSearchRangePPM));
          // multi-local-max
          //start
          if (rightPos == InvalidVal){
            if (count >= MIN_CHECK) {
              localpepMax[col_for_max].add(pos_for_max);
              if (prev_value_int_max!=InvalidVal && pit_intval < value_int_max * PercentageOfMax) {
                localpepMax[prev_col_for_max].add(prev_pos_for_max);
              }
            }
            count = 1;
            break;
          }
          else if (intMap[rightCol][rightPos] < prev_val) {
            if (value_int_max <= prev_value_int_max) {
              value_int_max = prev_value_int_max;
              col_for_max = prev_col_for_max;
              pos_for_max = prev_pos_for_max;
            }
            else if (prev_value_int_max!=InvalidVal && pit_intval < value_int_max * PercentageOfMax) { // the second max, value_int_max > prev_value_int_max
              localpepMax[prev_col_for_max].add(prev_pos_for_max);
            }
            prev_value_int_max = InvalidVal; // turned off when not necessary
            isDecrease = true;
          }
          else if (isDecrease) { // cur_val >= prev_val
            if (prev_val < value_int_max * PercentageOfMax) { // prev_val is lowest intensity point
              localpepMax[col_for_max].add(pos_for_max);
            } else {
              prev_col_for_max = col_for_max;
              prev_pos_for_max = pos_for_max;
              prev_value_int_max = value_int_max;
            }
            value_int_max = intMap[rightCol][rightPos];
            col_for_max = rightCol;
            pos_for_max = rightPos;
            pit_intval = prev_val;
            isDecrease = false;
          }
          if (intMap[rightCol][rightPos] > value_int_max) {
            value_int_max = intMap[rightCol][rightPos];
            col_for_max = rightCol;
            pos_for_max = rightPos;
          }
          // update
          last_col = rightCol;
          last_pos = rightPos;
          prev_val = intMap[rightCol][rightPos];
          count++;
          //end

          // turn the visited to invalid
          intMap[rightCol][rightPos] = InvalidVal;

        }
        // add last local max position;
        if (rightPos == IDX && count >= MIN_CHECK) {
          localpepMax[col_for_max].add(pos_for_max);
        }
      }
    }
    // sort
    for (int i = 0; i < IDX; i++) {
      Collections.sort(localpepMax[i]);
    }
  }


  /**
   * Search for isotopes of all charge states based on local maximal intensities; calculate
   * intensity shape function to measure feature quality, record only when isotopes >= 2
   *
   * @param window_size the max distance to search isotopes
   * @param num_of_z z values: 1 ... z
   */
  private void searchIsotopeWithSkipScan(int window_size, int num_of_z) {
    scoreIntShape = new ArrayList<>();
    scorePepAll = new ArrayList<>();
    scorePepAllPos = new ArrayList<>();
    scoreZVal = new ArrayList<>();
    intCluster = new ArrayList<>(); // Cluster of intensities corresponding to feature cluster
    intSum = new ArrayList<>();
    scoreIsoNum = new ArrayList<>();
    intensityPercentage = new ArrayList<>();
    int groupID = 0; // Corresponding to scorePepAll index

    System.out.println("IDX:"+IDX);
    System.out.println("mzMap:"+mzMap.length);
    // Search for feature all z
    for (int z = num_of_z; z >= 1; z--) {//z from max to 1
      List<Integer>[] maxMap = new List[IDX];//IDX is number of scan
      for (int i = 0; i < IDX; i++) {
        List<Integer> copy = new ArrayList<>(localpepMax[i]);
        maxMap[i] = copy;//get the local maximal intensity
      }
//      for (int i = 0; i < IDX - window_size + 1; i++) {//window_size is the max distance to search isotopes,
      for (int i = 0; i < IDX - window_size*iSkipScanMS1 + 1; i++) {//window_size is the max distance to search isotopes,
        System.out.println("i:"+i+" window_size:"+window_size+"  iskipscanMS1:"+iSkipScanMS1+" IDX:"+IDX);


        for (int k = 0; k < maxMap[i].size(); k++) {
          if (maxMap[i].get(k) == InvalidVal) {
            continue;
          }
          //when the intensity has a value
          double cur_mz = mzMap[i][maxMap[i].get(k)]; // M/Z value of current point
          double center_mz = cur_mz;
          List<Pair<Integer, Integer>> tempPos =
                  new ArrayList<>(); // All positions found in terms of localPepMax. {i, k}
          tempPos.add(new Pair<>(i, k));

          // Intensities of current peak and peaks within window size
          List<Double> curInt = new ArrayList<>();
          List<Double> maxInt = new ArrayList<>(); // max intensity for an isotope

          // Intensities for a whole cluster in increasing order of mz
          List<List<Double>> clusterInt = new ArrayList<>(); // intensities for a whole cluster
          curInt.add(intensityMap[i][maxMap[i].get(k)]);
          maxInt.add(intensityMap[i][maxMap[i].get(k)]);

          // Looking for intensity
          // Recorded as 0, if no reads
          // check next windows_size scan, whether there is mz same as curmz
//          for (int p = 1; p < window_size; p++) {
          for (int p = iSkipScanMS1; p < window_size*iSkipScanMS1; p+=iSkipScanMS1) {
            int pos_search =
                    searchRange(
                            mzMap[i + p], cur_mz * (1 - mzSearchRangePPM), cur_mz * (1 + mzSearchRangePPM));
            if (pos_search != InvalidVal) {
              curInt.add(intensityMap[i + p][pos_search]);
            } else {
              curInt.add(0.0);
            }
          }
          clusterInt.add(curInt);

          //look for isotope with z
          double lowrange = cur_mz * (1 - ppm) - (C13 - C12) / z;
          double hirange = lowrange + 2 * cur_mz * ppm;
          int lowerPos = InvalidVal;
          while (true) {
            int j = iSkipScanMS1;
            List<Double> newInt = new ArrayList<>();
//            for (; j < window_size; j++) {
            for (; j < window_size*iSkipScanMS1; j+=iSkipScanMS1) {
              lowerPos = searchPos(maxMap[i + j], mzMap[i + j], lowrange, hirange);
              if (lowerPos != InvalidVal) {
                cur_mz = mzMap[i + j][maxMap[i + j].get(lowerPos)];
                lowrange = cur_mz * (1 - ppm) - (C13 - C12) / z;
                hirange = lowrange + 2 * cur_mz * ppm;
                tempPos.add(0, new Pair<>(i + j, lowerPos));
                break;
              }
            }
            // If no point found below, break
            // If found, record intensity for the feature
            if (lowerPos == InvalidVal) {
              break;
            } else {
              for (int u = 0; u < window_size*iSkipScanMS1; u+=iSkipScanMS1) {
//              for (int u = 0; u < window_size; u++) {
                if (u == j) {
                  maxInt.add(0, intensityMap[i + j][maxMap[i + j].get(lowerPos)]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(lowerPos)]);
                  continue;
                }
                int pos_search =
                        searchRange(
                                mzMap[i+u], cur_mz * (1-mzSearchRangePPM), cur_mz * (1+mzSearchRangePPM));
                if (pos_search != -1) {
                  newInt.add(intensityMap[i+u][pos_search]);
                } else  {
                  newInt.add(0.00);
                }
              }
            }
            clusterInt.add(0, newInt);
          }

          // Check above the center peptide for isotopes
          lowrange = center_mz * (1 - ppm) + (C13 - C12) / z;
          hirange = lowrange + 2 * center_mz * ppm;
          int higherPos = InvalidVal;
          while (true) {
            int j = 0;
            List<Double> newInt = new ArrayList<>();
            for (; j < window_size*iSkipScanMS1; j+=iSkipScanMS1) {
//            for (; j < window_size; j++) {
              higherPos = searchPos(maxMap[i + j], mzMap[i + j], lowrange, hirange);
              if (higherPos != InvalidVal) {
                cur_mz = mzMap[i + j][maxMap[i + j].get(higherPos)];
                lowrange = cur_mz * (1 - ppm) + (C13 - C12) / z;
                hirange = lowrange + 2 * center_mz * ppm;
                tempPos.add(new Pair<>(i + j, higherPos));
                break;
              }
            }
            // If no point found above, break
            // If found, record intensity for the feature
            if (higherPos == InvalidVal) {
              break;
            } else {
//              for (int u = 0; u < window_size; u++) {
              for (int u = 0; u < window_size*iSkipScanMS1; u+=iSkipScanMS1) {
                if (u == j) {
                  maxInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos)]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos)]);
                  continue;
                }
                int pos_search =
                        searchRange(
                                mzMap[i + u], cur_mz * (1 - mzSearchRangePPM), cur_mz * (1 + mzSearchRangePPM));
                if (pos_search != InvalidVal) {
                  newInt.add(intensityMap[i + u][pos_search]);
                } else {
                  newInt.add(0.0);
                }
              }
            }
            clusterInt.add(newInt);
          }
          int num_pep = tempPos.size();

          // Check if there are >=MIN_CHECK isotopes and intensity of major isotope is >=intensityThreshold
          if (num_pep >= MIN_CHECK && maxInt.get(0) >= intensityThreshold) {
            double intShape = 0;

            // Calculate intensity shape score
            for (int t = 0; t < tempPos.size() - 1; t++) {
//              System.out.println(clusterInt.get(t).size()+": "+clusterInt.get(t + 1).size());
              intShape += addScore(clusterInt.get(t), clusterInt.get(t + 1));
            }
            intShape /= tempPos.size() - 1;
            // Calculate sum of intensity for all lines
            List<Double> intSumLine = new ArrayList<>();
            for (int t = 0; t < clusterInt.size(); t++) {
              double intSum = 0;
              for (int v = 0; v < clusterInt.get(t).size(); v++) {
                intSum += clusterInt.get(t).get(v);
              }
              intSumLine.add(intSum);
            }
            // Turn local max into invalid
            // Turn tempPos into real position
            List<Pair<Double, Double>> realPos = new ArrayList<>();
            List<Pair<Integer, Integer>> virtualPos = new ArrayList<>();
            for (int l = 0; l < tempPos.size(); l++) {
              // Turn local max into invalid
              maxMap[tempPos.get(l).getL()].set(tempPos.get(l).getR(), InvalidVal);
              // Turn tempPos into real position
              realPos.add(
                      new Pair<>(
                              RT[tempPos.get(l).getL()],
                              mzMap[tempPos.get(l).getL()][
                                      localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR())]));
              virtualPos.add(
                      new Pair<>(
                              tempPos.get(l).getL(),
                              localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR())));

            }

            // calculate intensity percentage of local area
            double lowMZ = realPos.get(0).getR() - IntensityBound;
            double hiMZ = realPos.get(realPos.size()-1).getR() + IntensityBound;
            double areaIntSum = 0;
            double isoIntSum = 0;
//            for (int m = 0; m < window_size; m++) {
            for (int m = 0; m < window_size*iSkipScanMS1; m+=iSkipScanMS1) {
              areaIntSum += intensityColSum(mzMap[i+m], intensityMap[i+m], lowMZ, hiMZ);
            }
            for (int m = 0; m < intSumLine.size(); m++) {
              isoIntSum += intSumLine.get(m);
            }
            double isoIntPercentage = isoIntSum/areaIntSum;

            groupID++;
            scorePepAll.add(realPos);
            scorePepAllPos.add(virtualPos);
            scoreIntShape.add(intShape);
            scoreZVal.add(z);
            scoreIsoNum.add(tempPos.size());
            intCluster.add(intSumLine);
            intSum.add(isoIntSum/window_size);
            intensityPercentage.add(isoIntPercentage);
          }
        }
      }
    }
    System.out.println("All peptide found before deleting duplicates: " + groupID);
  }

  /**
   * Search for isotopes of all charge states based on local maximal intensities; calculate
   * intensity shape function to measure feature quality, record only when isotopes >= 2
   *
   * @param window_size the max distance to search isotopes
   * @param num_of_z z values: 1 ... z
   */
  private void searchIsotope(int window_size, int num_of_z) {
    scoreIntShape = new ArrayList<>();
    scorePepAll = new ArrayList<>();
    scorePepAllPos = new ArrayList<>();
    scoreZVal = new ArrayList<>();
    intCluster = new ArrayList<>(); // Cluster of intensities corresponding to feature cluster
    intSum = new ArrayList<>();
    scoreIsoNum = new ArrayList<>();
    intensityPercentage = new ArrayList<>();
    int groupID = 0; // Corresponding to scorePepAll index

    // Search for feature all z
    for (int z = num_of_z; z >= 1; z--) {//z from max to 1
      List<Integer>[] maxMap = new List[IDX];//IDX is number of scan
      for (int i = 0; i < IDX; i++) {
        List<Integer> copy = new ArrayList<>(localpepMax[i]);
        maxMap[i] = copy;//get the local maximal intensity
      }
      for (int i = 0; i < IDX - window_size + 1; i++) {//window_size is the max distance to search isotopes,
        for (int k = 0; k < maxMap[i].size(); k++) {
          if (maxMap[i].get(k) == InvalidVal) {
            continue;
          }
          //when the intensity has a value
          double cur_mz = mzMap[i][maxMap[i].get(k)]; // M/Z value of current point
          double center_mz = cur_mz;
          List<Pair<Integer, Integer>> tempPos =
                  new ArrayList<>(); // All positions found in terms of localPepMax. {i, k}
          tempPos.add(new Pair<>(i, k));

          // Intensities of current peak and peaks within window size
          List<Double> curInt = new ArrayList<>();
          List<Double> maxInt = new ArrayList<>(); // max intensity for an isotope

          // Intensities for a whole cluster in increasing order of mz
          List<List<Double>> clusterInt = new ArrayList<>(); // intensities for a whole cluster
          curInt.add(intensityMap[i][maxMap[i].get(k)]);
          maxInt.add(intensityMap[i][maxMap[i].get(k)]);

          // Looking for intensity
          // Recorded as 0, if no reads
          // check next windows_size scan, whether there is mz same as curmz
          for (int p = 1; p < window_size; p++) {
            int pos_search =
                    searchRange(
                            mzMap[i + p], cur_mz * (1 - mzSearchRangePPM), cur_mz * (1 + mzSearchRangePPM));
            if (pos_search != InvalidVal) {
              curInt.add(intensityMap[i + p][pos_search]);
            } else {
              curInt.add(0.0);
            }
          }
          clusterInt.add(curInt);

          //look for isotope with z
          double lowrange = cur_mz * (1 - ppm) - (C13 - C12) / z;
          double hirange = lowrange + 2 * cur_mz * ppm;
          int lowerPos = InvalidVal;
          while (true) {
            int j = 1;
            List<Double> newInt = new ArrayList<>();
            for (; j < window_size; j++) {
              lowerPos = searchPos(maxMap[i + j], mzMap[i + j], lowrange, hirange);
              if (lowerPos != InvalidVal) {
                cur_mz = mzMap[i + j][maxMap[i + j].get(lowerPos)];
                lowrange = cur_mz * (1 - ppm) - (C13 - C12) / z;
                hirange = lowrange + 2 * cur_mz * ppm;
                tempPos.add(0, new Pair<>(i + j, lowerPos));
                break;
              }
            }
            // If no point found below, break
            // If found, record intensity for the feature
            if (lowerPos == InvalidVal) {
              break;
            } else {
              for (int u = 0; u < window_size; u++) {
                if (u == j) {
                  maxInt.add(0, intensityMap[i + j][maxMap[i + j].get(lowerPos)]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(lowerPos)]);
                  continue;
                }
                int pos_search =
                        searchRange(
                                mzMap[i+u], cur_mz * (1-mzSearchRangePPM), cur_mz * (1+mzSearchRangePPM));
                if (pos_search != -1) {
                  newInt.add(intensityMap[i+u][pos_search]);
                } else  {
                  newInt.add(0.00);
                }
              }
            }
            clusterInt.add(0, newInt);
          }

          // Check above the center peptide for isotopes
          lowrange = center_mz * (1 - ppm) + (C13 - C12) / z;
          hirange = lowrange + 2 * center_mz * ppm;
          int higherPos = InvalidVal;
          while (true) {
            int j = 0;
            List<Double> newInt = new ArrayList<>();
            for (; j < window_size; j++) {
              higherPos = searchPos(maxMap[i + j], mzMap[i + j], lowrange, hirange);
              if (higherPos != InvalidVal) {
                cur_mz = mzMap[i + j][maxMap[i + j].get(higherPos)];
                lowrange = cur_mz * (1 - ppm) + (C13 - C12) / z;
                hirange = lowrange + 2 * center_mz * ppm;
                tempPos.add(new Pair<>(i + j, higherPos));
                break;
              }
            }
            // If no point found above, break
            // If found, record intensity for the feature
            if (higherPos == InvalidVal) {
              break;
            } else {
              for (int u = 0; u < window_size; u++) {
                if (u == j) {
                  maxInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos)]);
                  newInt.add(intensityMap[i + j][maxMap[i + j].get(higherPos)]);
                  continue;
                }
                int pos_search =
                        searchRange(
                                mzMap[i + u], cur_mz * (1 - mzSearchRangePPM), cur_mz * (1 + mzSearchRangePPM));
                if (pos_search != InvalidVal) {
                  newInt.add(intensityMap[i + u][pos_search]);
                } else {
                  newInt.add(0.0);
                }
              }
            }
            clusterInt.add(newInt);
          }
          int num_pep = tempPos.size();

          // Check if there are >=MIN_CHECK isotopes and intensity of major isotope is >=intensityThreshold
          if (num_pep >= MIN_CHECK && maxInt.get(0) >= intensityThreshold) {
            double intShape = 0;

            // Calculate intensity shape score
            for (int t = 0; t < tempPos.size() - 1; t++) {
              intShape += addScore(clusterInt.get(t), clusterInt.get(t + 1));
            }
            intShape /= tempPos.size() - 1;
            // Calculate sum of intensity for all lines
            List<Double> intSumLine = new ArrayList<>();
            for (int t = 0; t < clusterInt.size(); t++) {
              double intSum = 0;
              for (int v = 0; v < clusterInt.get(t).size(); v++) {
                intSum += clusterInt.get(t).get(v);
              }
              intSumLine.add(intSum);
            }
            // Turn local max into invalid
            // Turn tempPos into real position
            List<Pair<Double, Double>> realPos = new ArrayList<>();
            List<Pair<Integer, Integer>> virtualPos = new ArrayList<>();
            for (int l = 0; l < tempPos.size(); l++) {
              // Turn local max into invalid
              maxMap[tempPos.get(l).getL()].set(tempPos.get(l).getR(), InvalidVal);
              // Turn tempPos into real position
              realPos.add(
                      new Pair<>(
                              RT[tempPos.get(l).getL()],
                              mzMap[tempPos.get(l).getL()][
                                      localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR())]));
              virtualPos.add(
                      new Pair<>(
                              tempPos.get(l).getL(),
                              localpepMax[tempPos.get(l).getL()].get(tempPos.get(l).getR())));

            }

            // calculate intensity percentage of local area
            double lowMZ = realPos.get(0).getR() - IntensityBound;
            double hiMZ = realPos.get(realPos.size()-1).getR() + IntensityBound;
            double areaIntSum = 0;
            double isoIntSum = 0;
            for (int m = 0; m < window_size; m++) {
              areaIntSum += intensityColSum(mzMap[i+m], intensityMap[i+m], lowMZ, hiMZ);
            }
            for (int m = 0; m < intSumLine.size(); m++) {
              isoIntSum += intSumLine.get(m);
            }
            double isoIntPercentage = isoIntSum/areaIntSum;

            groupID++;
            scorePepAll.add(realPos);
            scorePepAllPos.add(virtualPos);
            scoreIntShape.add(intShape);
            scoreZVal.add(z);
            scoreIsoNum.add(tempPos.size());
            intCluster.add(intSumLine);
            intSum.add(isoIntSum/window_size);
            intensityPercentage.add(isoIntPercentage);
          }
        }
      }
    }
    System.out.println("All peptide found before deleting duplicates: " + groupID);
  }

  /**
   * Search for isotopes of all charge states based on raw signal data before local maximal intensities;
   * calculate
   *
   * @param window_size the max distance to search isotopes
   * @param num_of_z z values: 1 ... z
   */
  //@author JCZhong
  public void searchIsotopeWithSignalWithSkipSpan(int window_size, int num_of_z) {
    //Check all raw mz value, find all possible isotope signal in raw data
    //for 1 charge, mz + 1 with ppm
    //for 2 charge, mz + 0.5 with ppm
    //for N charge, mz + 1/N with ppm
    // Search max near one m/z vertically[current mass value,+N(N=10)] default charge value is 10
//    System.out.println("idx = " + IDX);
    ProgressBar pb = new ProgressBar("Progress", IDX);
    pb.start();
    int id=1;
    String spectrumFolder = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/";

    FileWriter fileWriterMS1Isotope = null;
    try {
//      fileWriterMS1Isotope = new FileWriter(spectrumFolder + "ms1istopeRemoveHighChargeCover.csv");
      fileWriterMS1Isotope = new FileWriter(spectrumFolder + "ms1istopeRemoveHighChargeCoverR03.csv");
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    BufferedWriter bufferedWriterMSIsotope = new BufferedWriter(fileWriterMS1Isotope);
    try {
//      bufferedWriterMSIsotope.write("timeindex\ttime\tmz\tcharges\tistopes\n");
//      bufferedWriterMSIsotope.write("timeindex\ttime\tmz\tcharges\tistopes\n");

      String header =
              "id"
                      + '\t'
                      + "mz"
                      + '\t'
                      + "rt"
                      + '\t'
                      + "z"
                      + '\t'
                      + "isotope_num"
                      + '\t'
                      + "scan_num"
                      + '\t'
                      + "quantification_peaks_area"
                      + '\n';
      bufferedWriterMSIsotope.write(header);

//      bufferedWriterMSIsotope.write("timeindex\trt\tmz\tintensity\tz\tistopes\n");
    } catch (IOException e) {
      throw new RuntimeException(e);
    }

    List<SignalIsotope> liSigIsotope = new ArrayList<>();
    int num_of_Isotope = 8;//
    for (int width = 0; width < IDX; width++) {
      pb.step();
      int scanMassCount = mzMap[width].length;
      for (int height = 0; height < scanMassCount-1; height++) {
        SignalIsotope sIsotope = null;
        int iCurPos = height;
//        double[] dLowMzOfIsotopeWithZ = new double[num_of_z];
//        double[] dHighMzofIsotopeWithZ = new double[num_of_z];
//        for (int iZ = 1; iZ <= num_of_z; iZ++) {
//          dLowMzOfIsotopeWithZ[iZ-1] = (mzMap[width][height] + 1.0/iZ) * (1 - ppm);//1 charge with an isotope
//          dHighMzofIsotopeWithZ[iZ-1] = (mzMap[width][height] + 1.0 * num_of_Isotope/iZ) * (1 + ppm);//1 charge with number of isotope
//        }

//        for (int iZ = 1; iZ <= num_of_z; iZ++) {
//        if (width==1)
//          System.out.println("test");


//        if(mzMap[width][height]>316.156 && mzMap[width][height]<316.157)
//          System.out.println("test");

        for (int iZ = num_of_z; iZ >= 1; iZ--) {
          List<Integer> listZIsotopePos = null;
          List<Double> listZIsotopeValue = null;

          int iIsotpe = 1;
          boolean bzHighChargeCover = false;
          boolean bzFindNew = false;

          double dHighMzofIsotopeWithZ = mzMap[width][height] * (1 + ppm) + num_of_Isotope * (C13 - C12) / iZ;
//        double dHighMzofIsotopeWithZ = (mzMap[width][height] + 1.0 * num_of_Isotope) * (1 + ppm);//1 charge with number of isotope
          int iNextHighIsotopePos = getNearpos(mzMap[width], dHighMzofIsotopeWithZ, iCurPos, scanMassCount);//the maximum pos
          if (iNextHighIsotopePos + 1 <= scanMassCount) {
            iNextHighIsotopePos++;//the pos  greater than the dHighMzofIsotopeWithZ
          }
          int iCurrentNearPos = 0;
          int iCurrentZ_num_Isotope =0;
          for (; iIsotpe <= num_of_Isotope; iIsotpe++) {//find the nearest pos
            iCurrentNearPos = 0;
            bzHighChargeCover = false;
            //IF Z is half of num_of_z,we should j
            if (iZ <= num_of_z / 2) {
              int tmpZ = iZ;
              int tmpIsotope = iIsotpe;
              while (tmpZ+iZ <= num_of_z) {

                tmpZ += iZ;
                tmpIsotope = iIsotpe * tmpZ / iZ;

                if (sIsotope!= null && sIsotope.lliIstopePos!=null && sIsotope.lliIstopePos.get(tmpZ) != null) {
                  if (sIsotope.lliIstopePos.get(tmpZ).size() >= tmpIsotope) {
//                        bzAreadyFound = sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1) == iCurrentNearPos;
                    bzHighChargeCover = true;
                    if (listZIsotopePos == null) {
                      listZIsotopePos = new ArrayList<>();
                      listZIsotopeValue = new ArrayList<>();
                    }
                    listZIsotopePos.add(sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1));
                    listZIsotopeValue.add(mzMap[width][sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1)]);
                  }
                }
              }
            }
            if (!bzHighChargeCover) {//没有覆盖
              double dRealMzOfIsotopeWithZ = mzMap[width][height];
//            double dLowMzOfIsotopeWithZ = (mzMap[width][height] + iIsotpe * 1.0/ iZ );//1 charge with an isotope
              int iNextIsotopePos = getNearpos(mzMap[width], dRealMzOfIsotopeWithZ + iIsotpe * (C13 - C12) / iZ, iCurPos + 1, iNextHighIsotopePos);

//            if( (mzMap[width][iNextIsotopePos] >= dRealMzOfIsotopeWithZ*(1-ppm) + num_of_Isotope * (C13 - C12) / iZ) &&
//                    (mzMap[width][iNextIsotopePos] <= dRealMzOfIsotopeWithZ*(1+ppm) + num_of_Isotope * (C13 - C12) / iZ))
//            {
              iCurrentNearPos = isInPPM(mzMap[width], dRealMzOfIsotopeWithZ, iNextIsotopePos, iIsotpe, iZ);
            }
            if (iCurrentNearPos > 0) {

              if (sIsotope == null) {
                sIsotope = new SignalIsotope();
                sIsotope.iRTpos = width;
                sIsotope.iMzPos = height;
                sIsotope.mzvalue = mzMap[width][height];
//                  sIsotope.liCharge = new HashMap<>();
                sIsotope.lliIstopePos = new HashMap<>();
                sIsotope.lliIstopeValue = new HashMap<>();
              }
              if (listZIsotopePos == null) {
                listZIsotopePos = new ArrayList<>();
                listZIsotopeValue = new ArrayList<>();
              }
              listZIsotopePos.add(iCurrentNearPos);
              listZIsotopeValue.add(mzMap[width][iCurrentNearPos]);
              iCurrentZ_num_Isotope = iIsotpe;
              bzFindNew = true;

            } else if(!bzHighChargeCover) {
              break;
            }

//            }
//            else {
//              break;
//            }
          }
          if (bzFindNew)//NO COVERED BY HIGH CHARGE AND FIND NEW
          {


//            sIsotope.liCharge.add(iZ);
//            sIsotope.lliIstopePos.add(listZIsotopePos);
//            sIsotope.lliIstopeValue.add(listZIsotopeValue);
            sIsotope.lliIstopePos.put(iZ,listZIsotopePos);
            sIsotope.lliIstopeValue.put(iZ,listZIsotopeValue);

            try {

//              bufferedWriterMSIsotope.write(id++ +"\t"+mzMap[width][height]+"\t"+RT[width]+"\t"+iZ+
//                      "\t"+iCurrentZ_num_Isotope+"\t0\t0\t0\t"+width+"\t"+ intensityMap[width][height]+"\t"+ intensityMap[width][height]+
//                      "\t0\t0\t"+mzMap[width][height]+"\t"+RT[width]+"\t"+intensityMap[width][height]+"\n");
              bufferedWriterMSIsotope.write(id++ +"\t"+mzMap[width][height]+"\t"+RT[width]+"\t"+iZ+
                      "\t"+iCurrentZ_num_Isotope+"\t"+width+"\t"+ intensityMap[width][height]+"\n");
            } catch (IOException e) {
              throw new RuntimeException(e);
            }
          }

        }

        if (sIsotope!=null)
        {
          liSigIsotope.add(sIsotope);
        }
      }
    }
    try {
      bufferedWriterMSIsotope.flush();
      bufferedWriterMSIsotope.close();;

    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    pb.stop();
//    System.out.println("end");
  }
  /**
   * Search for isotopes of all charge states based on raw signal data before local maximal intensities;
   * calculate
   *
   * @param window_size the max distance to search isotopes
   * @param num_of_z z values: 1 ... z
   */
  //@author JCZhong
  public void searchIsotopeWithSignal(int window_size, int num_of_z) {
    //Check all raw mz value, find all possible isotope signal in raw data
    //for 1 charge, mz + 1 with ppm
    //for 2 charge, mz + 0.5 with ppm
    //for N charge, mz + 1/N with ppm
    // Search max near one m/z vertically[current mass value,+N(N=10)] default charge value is 10
//    System.out.println("idx = " + IDX);
    ProgressBar pb = new ProgressBar("Progress", IDX);
    pb.start();
    int id=1;
    String spectrumFolder = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/mouse/";

    FileWriter fileWriterMS1Isotope = null;
    try {
//      fileWriterMS1Isotope = new FileWriter(spectrumFolder + "ms1istopeRemoveHighChargeCover.csv");
      fileWriterMS1Isotope = new FileWriter(spectrumFolder + "ms1istopeRemoveHighChargeCoverMouseR02.csv");
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    BufferedWriter bufferedWriterMSIsotope = new BufferedWriter(fileWriterMS1Isotope);
    try {
//      bufferedWriterMSIsotope.write("timeindex\ttime\tmz\tcharges\tistopes\n");
//      bufferedWriterMSIsotope.write("timeindex\ttime\tmz\tcharges\tistopes\n");

      String header =
              "id"
                      + '\t'
                      + "mz"
                      + '\t'
                      + "rt"
                      + '\t'
                      + "z"
                      + '\t'
                      + "isotope_num"
                      + '\t'
                      + "scan_num"
                      + '\t'
                      + "quantification_peaks_area"
                      + '\n';
      bufferedWriterMSIsotope.write(header);

//      bufferedWriterMSIsotope.write("timeindex\trt\tmz\tintensity\tz\tistopes\n");
    } catch (IOException e) {
      throw new RuntimeException(e);
    }

    List<SignalIsotope> liSigIsotope = new ArrayList<>();
    int num_of_Isotope = 8;//
    for (int width = 0; width < IDX; width++) {
      pb.step();
      int scanMassCount = mzMap[width].length;
      for (int height = 0; height < scanMassCount-1; height++) {
        SignalIsotope sIsotope = null;
        int iCurPos = height;
//        double[] dLowMzOfIsotopeWithZ = new double[num_of_z];
//        double[] dHighMzofIsotopeWithZ = new double[num_of_z];
//        for (int iZ = 1; iZ <= num_of_z; iZ++) {
//          dLowMzOfIsotopeWithZ[iZ-1] = (mzMap[width][height] + 1.0/iZ) * (1 - ppm);//1 charge with an isotope
//          dHighMzofIsotopeWithZ[iZ-1] = (mzMap[width][height] + 1.0 * num_of_Isotope/iZ) * (1 + ppm);//1 charge with number of isotope
//        }

//        for (int iZ = 1; iZ <= num_of_z; iZ++) {
//        if (width==1)
//          System.out.println("test");


//        if(mzMap[width][height]>316.156 && mzMap[width][height]<316.157)
//          System.out.println("test");

        for (int iZ = num_of_z; iZ >= 1; iZ--) {
          List<Integer> listZIsotopePos = null;
          List<Double> listZIsotopeValue = null;

          int iIsotpe = 1;
          boolean bzHighChargeCover = false;
          boolean bzFindNew = false;

          double dHighMzofIsotopeWithZ = mzMap[width][height] * (1 + ppm) + num_of_Isotope * (C13 - C12) / iZ;
//        double dHighMzofIsotopeWithZ = (mzMap[width][height] + 1.0 * num_of_Isotope) * (1 + ppm);//1 charge with number of isotope
          int iNextHighIsotopePos = getNearpos(mzMap[width], dHighMzofIsotopeWithZ, iCurPos, scanMassCount);//the maximum pos
          if (iNextHighIsotopePos + 1 <= scanMassCount) {
            iNextHighIsotopePos++;//the pos  greater than the dHighMzofIsotopeWithZ
          }

          int iCurrentNearPos = 0;
          int iCurrentZ_num_Isotope =0;
          for (; iIsotpe <= num_of_Isotope; iIsotpe++) {//find the nearest pos
            iCurrentNearPos = 0;
            bzHighChargeCover = false;
            //IF Z is half of num_of_z,we should j
            if (iZ <= num_of_z / 2) {
              int tmpZ = iZ;
              int tmpIsotope = iIsotpe;
              while (tmpZ+iZ <= num_of_z) {

                tmpZ += iZ;
                tmpIsotope = iIsotpe * tmpZ / iZ;

                if (sIsotope!= null && sIsotope.lliIstopePos!=null && sIsotope.lliIstopePos.get(tmpZ) != null) {
                  if (sIsotope.lliIstopePos.get(tmpZ).size() >= tmpIsotope) {
//                        bzAreadyFound = sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1) == iCurrentNearPos;
                    bzHighChargeCover = true;
                    if (listZIsotopePos == null) {
                      listZIsotopePos = new ArrayList<>();
                      listZIsotopeValue = new ArrayList<>();
                    }
                    listZIsotopePos.add(sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1));
                    listZIsotopeValue.add(mzMap[width][sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1)]);
                  }
                }
              }
            }
            if (!bzHighChargeCover) {//没有覆盖
              double dRealMzOfIsotopeWithZ = mzMap[width][height];
//            double dLowMzOfIsotopeWithZ = (mzMap[width][height] + iIsotpe * 1.0/ iZ );//1 charge with an isotope
              int iNextIsotopePos = getNearpos(mzMap[width], dRealMzOfIsotopeWithZ + iIsotpe * (C13 - C12) / iZ, iCurPos + 1, iNextHighIsotopePos);

//            if( (mzMap[width][iNextIsotopePos] >= dRealMzOfIsotopeWithZ*(1-ppm) + num_of_Isotope * (C13 - C12) / iZ) &&
//                    (mzMap[width][iNextIsotopePos] <= dRealMzOfIsotopeWithZ*(1+ppm) + num_of_Isotope * (C13 - C12) / iZ))
//            {
              iCurrentNearPos = isInPPM(mzMap[width], dRealMzOfIsotopeWithZ, iNextIsotopePos, iIsotpe, iZ);
            }
            if (iCurrentNearPos > 0) {

              if (sIsotope == null) {
                sIsotope = new SignalIsotope();
                sIsotope.iRTpos = width;
                sIsotope.iMzPos = height;
                sIsotope.mzvalue = mzMap[width][height];
//                  sIsotope.liCharge = new HashMap<>();
                sIsotope.lliIstopePos = new HashMap<>();
                sIsotope.lliIstopeValue = new HashMap<>();
              }
              if (listZIsotopePos == null) {
                listZIsotopePos = new ArrayList<>();
                listZIsotopeValue = new ArrayList<>();
              }
              listZIsotopePos.add(iCurrentNearPos);
              listZIsotopeValue.add(mzMap[width][iCurrentNearPos]);
              iCurrentZ_num_Isotope = iIsotpe;
              bzFindNew = true;

            } else if(!bzHighChargeCover) {
              break;
            }


          }
          if (bzFindNew)//NO COVERED BY HIGH CHARGE AND FIND NEW
          {



            sIsotope.lliIstopePos.put(iZ,listZIsotopePos);
            sIsotope.lliIstopeValue.put(iZ,listZIsotopeValue);

            try {


              bufferedWriterMSIsotope.write(id++ +"\t"+mzMap[width][height]+"\t"+RT[width]+"\t"+iZ+
                      "\t"+iCurrentZ_num_Isotope+"\t"+width+"\t"+ intensityMap[width][height]+"\n");
            } catch (IOException e) {
              throw new RuntimeException(e);
            }
          }

        }

        if (sIsotope!=null)
        {
          liSigIsotope.add(sIsotope);

        }


      }
    }
    try {
      bufferedWriterMSIsotope.flush();
      bufferedWriterMSIsotope.close();;

    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    pb.stop();
//    System.out.println("end");
  }
  /*
  if find then return the pos
  else return 0
   */
  private int isInPPM(double[] mzs, double dLowMzOfIsotopeWithZ, int iNextHighIsotopePos) {
    int  iInPPMNearestPos = 0;
    double maxMZWithPPM = dLowMzOfIsotopeWithZ*(1+ppm);
    double minMzWithPPM = dLowMzOfIsotopeWithZ*(1-ppm);
    double mindistance = minMzWithPPM - dLowMzOfIsotopeWithZ;
    for (;iNextHighIsotopePos<mzs.length && mzs[iNextHighIsotopePos]<=maxMZWithPPM;iNextHighIsotopePos++)
    {
      if(abs(mzs[iNextHighIsotopePos]-dLowMzOfIsotopeWithZ)<=mindistance) {
        mindistance = mzs[iNextHighIsotopePos] - dLowMzOfIsotopeWithZ;
        iInPPMNearestPos = iNextHighIsotopePos;
      }else if(iInPPMNearestPos>0)
      {
        break;
      }

    }
    return iInPPMNearestPos;

  }

  private int isInPPM(double[] mzs, double drealMzOfIsotopeWithZ, int iNextHighIsotopePos,int isotope, int z) {
    int  iInPPMNearestPos = 0;
    double maxMZWithPPM = drealMzOfIsotopeWithZ*(1+ppm) + isotope * (C13 - C12) / z;
    double minMzWithPPM = drealMzOfIsotopeWithZ*(1-ppm) + isotope * (C13 - C12) / z;
    double realValue = drealMzOfIsotopeWithZ + isotope * (C13 - C12) / z;

    double mindistance = abs(minMzWithPPM - realValue);
    for (;iNextHighIsotopePos<mzs.length && mzs[iNextHighIsotopePos]<=maxMZWithPPM;iNextHighIsotopePos++)
    {
      if(abs(mzs[iNextHighIsotopePos]-realValue)<=mindistance) {
        mindistance = mzs[iNextHighIsotopePos] - drealMzOfIsotopeWithZ;
        iInPPMNearestPos = iNextHighIsotopePos;
      }else if(iInPPMNearestPos>0)
      {
        break;
      }

    }
    return iInPPMNearestPos;

  }

  /*
  find the nearest pos(less than value and the nearest)
   */
  public static int getiPosBeginOfNearest(int length, int iPosBegin) {
    if (iPosBegin < 0) {
      int indexOfNearest = ~iPosBegin;

      if (indexOfNearest == length) {
        //from time is larger than all elements
        iPosBegin = indexOfNearest;
      } else if (indexOfNearest == 0) {
        // from time is less than first item
        iPosBegin = 0;
      } else {
        // from time is between (indexOfNearest - 1) and indexOfNearest
        iPosBegin = indexOfNearest;
      }
    }
    return iPosBegin;
  }
  /*
    find the nearest pos based on binary search
     */
  private int getNearpos(double[] ms2, double dMz,int istartpos ,int iendpos) {
    int iposFront = binarySearch(ms2, istartpos, iendpos/*ms2.length*/, dMz);//全新查找
    int inearpos = getiPosBeginOfNearest(ms2.length, iposFront);

    if (inearpos == ms2.length) {
      inearpos = inearpos - 1;
    } else if ((inearpos > 0) && (abs(ms2[inearpos - 1] - dMz) < abs(ms2[inearpos] - dMz))) {
      inearpos = inearpos - 1;
    }
    return inearpos;
  }

  private void isoDistributionScore() {
    scoreIsoDistr = new ArrayList<>();
    for (int cur_pep = 0; cur_pep < intCluster.size(); cur_pep++) {
      // Calculate vector score of experimental and theoretical isotope intensity
      float m_over_z = scorePepAll.get(cur_pep).get(0).getR().floatValue();
      float m = m_over_z * scoreZVal.get(cur_pep);
      int maxMass = (int) (m + (C13 - C12) * scoreIsoNum.get(cur_pep) / scoreZVal.get(cur_pep) + 1);
      TheoreticalIsotope iso =
              new TheoreticalIsotope(scoreIsoNum.get(cur_pep), maxMass, MassStep, MinRelHeight);
      float[] iso_shape = iso.computeShape(m);
      int IsoNum = iso_shape.length;
      double xi_sqr_sum = 0;
      double yi_sqr_sum = 0;
      double[] xi = new double[IsoNum];

      // Calculate isotope distribution score
      for (int i = 0; i < IsoNum; i++) {
        yi_sqr_sum += Math.pow(iso_shape[i], 2);
        xi[i] = intCluster.get(cur_pep).get(i);
        xi_sqr_sum += Math.pow(xi[i], 2);
      }

      double xi_yi_dot_sum = 0;
      for (int i = 0; i < IsoNum; i++) {
        xi_yi_dot_sum += xi[i] * iso_shape[i];
      }
      double iso_score = xi_yi_dot_sum / (Math.sqrt(xi_sqr_sum) * Math.sqrt(yi_sqr_sum));
      scoreIsoDistr.add(iso_score);
    }
  }
  private void detectTrailofPrecursorWithSkipScan(){
    trailsAll = new ArrayList<>();

    for (List<Pair<Integer, Integer>> cluster: scorePepAllPos) {
      List<XIC> trails = new ArrayList<>();

      for (Pair<Integer, Integer> isotope: cluster) {
        ArrayList<Double>
                ints = new ArrayList<>(),
                mzs = new ArrayList<>(),
                rts = new ArrayList<>();

        ints.add(intensityMap[isotope.getL()][isotope.getR()]);
        mzs.add(mzMap[isotope.getL()][isotope.getR()]);
        rts.add(RT[isotope.getL()]);

        double lowMZ = mzMap[isotope.getL()][isotope.getR()] * (1-mzSearchRangePPM);
        double hiMZ = mzMap[isotope.getL()][isotope.getR()] * (1+mzSearchRangePPM);
        int leftScan = isotope.getL() - iSkipScanMS1;
//        int leftScan = isotope.getL() - 1;
        int rightScan = isotope.getL() + iSkipScanMS1;
//        int rightScan = isotope.getL() + 1;
        double leftMaxIntensity = intensityMap[isotope.getL()][isotope.getR()];
        double rightMaxIntensity = intensityMap[isotope.getL()][isotope.getR()];
//        for (;leftScan >= 0 ; leftScan--){
        for (;leftScan >= 0 ; leftScan-=iSkipScanMS1){
          int leftPos = searchRange(mzMap[leftScan], lowMZ, hiMZ);
          if (leftPos == InvalidVal || intensityMap[leftScan][leftPos] > leftMaxIntensity) {
            break;
          }
          lowMZ = mzMap[leftScan][leftPos] * (1-mzSearchRangePPM);
          hiMZ = mzMap[leftScan][leftPos] * (1+mzSearchRangePPM);
          leftMaxIntensity = intensityMap[leftScan][leftPos];

          ints.add(0, intensityMap[leftScan][leftPos]);
          mzs.add(0, mzMap[leftScan][leftPos]);
          rts.add(0, RT[leftScan]);

        }
//        for (;rightScan < IDX ; rightScan++){
        for (;rightScan < IDX ; rightScan+=iSkipScanMS1){
          int rightPos = searchRange(mzMap[rightScan], lowMZ, hiMZ);
          if (rightPos == InvalidVal || intensityMap[rightScan][rightPos] > rightMaxIntensity) {
            break;
          }
          lowMZ = mzMap[rightScan][rightPos] * (1-mzSearchRangePPM);
          hiMZ = mzMap[rightScan][rightPos] * (1+mzSearchRangePPM);
          rightMaxIntensity = intensityMap[rightScan][rightPos];

          ints.add(intensityMap[rightScan][rightPos]);
          mzs.add(mzMap[rightScan][rightPos]);
          rts.add( RT[rightScan]);
        }
        if (leftScan < 0) { leftScan = 0; }
        if (rightScan == IDX) { rightScan = IDX - 1; }

        XIC trail = new XIC(ints, mzs, rts);
        trails.add(trail);
      }
      trailsAll.add(trails);
    }
  }
  private void detectTrailofPrecursor(){
    trailsAll = new ArrayList<>();

    for (List<Pair<Integer, Integer>> cluster: scorePepAllPos) {
      List<XIC> trails = new ArrayList<>();

      for (Pair<Integer, Integer> isotope: cluster) {
        ArrayList<Double>
                ints = new ArrayList<>(),
                mzs = new ArrayList<>(),
                rts = new ArrayList<>();

        ints.add(intensityMap[isotope.getL()][isotope.getR()]);
        mzs.add(mzMap[isotope.getL()][isotope.getR()]);
        rts.add(RT[isotope.getL()]);

        double lowMZ = mzMap[isotope.getL()][isotope.getR()] * (1-mzSearchRangePPM);
        double hiMZ = mzMap[isotope.getL()][isotope.getR()] * (1+mzSearchRangePPM);
        int leftScan = isotope.getL() - 1;
        int rightScan = isotope.getL() + 1;
        double leftMaxIntensity = intensityMap[isotope.getL()][isotope.getR()];
        double rightMaxIntensity = intensityMap[isotope.getL()][isotope.getR()];
        for (;leftScan >= 0 ; leftScan--){
          int leftPos = searchRange(mzMap[leftScan], lowMZ, hiMZ);
          if (leftPos == InvalidVal || intensityMap[leftScan][leftPos] > leftMaxIntensity) {
            break;
          }
          lowMZ = mzMap[leftScan][leftPos] * (1-mzSearchRangePPM);
          hiMZ = mzMap[leftScan][leftPos] * (1+mzSearchRangePPM);
          leftMaxIntensity = intensityMap[leftScan][leftPos];

          ints.add(0, intensityMap[leftScan][leftPos]);
          mzs.add(0, mzMap[leftScan][leftPos]);
          rts.add(0, RT[leftScan]);

        }
        for (;rightScan < IDX ; rightScan++){
          int rightPos = searchRange(mzMap[rightScan], lowMZ, hiMZ);
          if (rightPos == InvalidVal || intensityMap[rightScan][rightPos] > rightMaxIntensity) {
            break;
          }
          lowMZ = mzMap[rightScan][rightPos] * (1-mzSearchRangePPM);
          hiMZ = mzMap[rightScan][rightPos] * (1+mzSearchRangePPM);
          rightMaxIntensity = intensityMap[rightScan][rightPos];

          ints.add(intensityMap[rightScan][rightPos]);
          mzs.add(mzMap[rightScan][rightPos]);
          rts.add( RT[rightScan]);
        }
        if (leftScan < 0) { leftScan = 0; }
        if (rightScan == IDX) { rightScan = IDX - 1; }

        XIC trail = new XIC(ints, mzs, rts);
        trails.add(trail);
      }
      trailsAll.add(trails);
    }
  }

  public int searchRange(double[] array, double lowrange, double highrange) {
    int pos = InvalidVal;
    for (int i = 0; i < array.length; i++) {
      if (array[i] >= lowrange && array[i] <= highrange) {
        pos = i;
        break;
      }
      if (array[i] > highrange) {
        break;
      }
    }
    return pos;
  }
  // Given array of localmax positions, check in range
  // Return position at posArray
  public int searchPos(
          List<Integer> posArray, double[] mzArray, double lowrange, double highrange) {
    int pos = InvalidVal;
    for (int i = 0; i < posArray.size(); i++) {
      if (posArray.get(i) == InvalidVal) {
        continue;
      }
      if (mzArray[posArray.get(i)] >= lowrange && mzArray[posArray.get(i)] <= highrange) {
        pos = i;
        break;
      }
      if (mzArray[posArray.get(i)] > highrange) {
        break;
      }
    }
    return pos;
  }

  double addScore(List<Double> a, List<Double> b) {
    double dot_sum = 0;
    double a_sqr_sum = 0;
    double b_sqr_sum = 0;
    for (int i = 0; i < a.size(); i++) {
      dot_sum += a.get(i) * b.get(i);
      a_sqr_sum += a.get(i) * a.get(i);
      b_sqr_sum += b.get(i) * b.get(i);
    }
    return dot_sum / (Math.sqrt(a_sqr_sum) * Math.sqrt(b_sqr_sum));
  }

  public double intensityColSum(double[] array, double[] intensity, double lowrange, double highrange){
    double intSum = 0;
    for (int i = 0; i < array.length; i++) {
      if (array[i] >= lowrange && array[i] <= highrange) {
        intSum += intensity[i];
      }
      if (array[i] > highrange) {
        break;
      }
    }
    return intSum;
  }

  private void sortInfoListParallel(List<double[]> infoList, int index) {
    mergeSortParallel(infoList, index, 0, infoList.size() - 1);
    // debug
//    for (int i = 0; i < infoList.size() - 1; i++) {
//      if (!(infoList.get(i)[index] <= infoList.get(i + 1)[index])) {
//        System.out.println("WARNING: InfoList Sort Failed.");
//        System.exit(-1);
//      }
//    }
  }
  private  void mergeSortParallel(List<double[]> infoList, int index, int l, int r) {
    if (r > l) {
      int m = l + (r - l)/2;
      mergeSortParallel(infoList, index, l, m);
      mergeSortParallel(infoList, index, m + 1, r);
      mergeParallel(infoList, index, l, m, r);
    }
  }
  private  void mergeParallel(List<double[]> infoList, int index, int l, int m, int r) {
    int l_size = m - l + 1, r_size = r - m;
    ArrayList<double[]> l_half = new ArrayList<>(l_size);
    ArrayList<double[]> r_half = new ArrayList<>(r_size);
    for (int i = 0; i < l_size; i++) {
      l_half.add(infoList.get(l + i));
    }
    for (int j = 0; j < r_size; j++) {
      r_half.add(infoList.get(m + j + 1));
    }

    int i = 0, j = 0, k = l;
    while (i < l_size && j < r_size) {
      if (l_half.get(i)[index] < r_half.get(j)[index]) {
        infoList.set(k, l_half.get(i));
        i++;
      } else {
        infoList.set(k, r_half.get(j));
        j++;
      }
      k++;
    }
    while (j < r_size) {
      infoList.set(k, r_half.get(j));
      j++;
      k++;
    }
    while (i < l_size) {
      infoList.set(k, l_half.get(i));
      i++;
      k++;
    }
  }
  private void sortInfoList(ArrayList<double[]> infoList, int index) {
    mergeSort(infoList, index, 0, infoList.size() - 1);
    // debug
//    for (int i = 0; i < infoList.size() - 1; i++) {
//      if (!(infoList.get(i)[index] <= infoList.get(i + 1)[index])) {
//        System.out.println("WARNING: InfoList Sort Failed.");
//        System.exit(-1);
//      }
//    }
  }
  private  void mergeSort(ArrayList<double[]> infoList, int index, int l, int r) {
    if (r > l) {
      int m = l + (r - l)/2;
      mergeSort(infoList, index, l, m);
      mergeSort(infoList, index, m + 1, r);
      merge(infoList, index, l, m, r);
    }
  }
  private  void merge(ArrayList<double[]> infoList, int index, int l, int m, int r) {
    int l_size = m - l + 1, r_size = r - m;
    ArrayList<double[]> l_half = new ArrayList<>(l_size);
    ArrayList<double[]> r_half = new ArrayList<>(r_size);
    for (int i = 0; i < l_size; i++) {
      l_half.add(infoList.get(l + i));
    }
    for (int j = 0; j < r_size; j++) {
      r_half.add(infoList.get(m + j + 1));
    }

    int i = 0, j = 0, k = l;
    while (i < l_size && j < r_size) {
      if (l_half.get(i)[index] < r_half.get(j)[index]) {
        infoList.set(k, l_half.get(i));
        i++;
      } else {
        infoList.set(k, r_half.get(j));
        j++;
      }
      k++;
    }
    while (j < r_size) {
      infoList.set(k, r_half.get(j));
      j++;
      k++;
    }
    while (i < l_size) {
      infoList.set(k, l_half.get(i));
      i++;
      k++;
    }
  }

  void writeFeatures(String filename) throws IOException {
    File file = new File(filename);
    try {
      // create FileWriter object with file as parameter
      FileWriter outputfile = new FileWriter(file);
      // PrintWriter
      PrintWriter printWriter = new PrintWriter(outputfile);
      // add header
      String header =
              "id"
                      + '\t'
                      + "mz"
                      + '\t'
                      + "rt"
                      + '\t'
                      + "z"
                      + '\t'
                      + "isotope_num"
                      + '\t'
                      + "intensity_shape_score"
                      + '\t'
                      + "isotope_distribution_score"
                      + '\t'
                      + "intensity_area_percentage"
                      + '\t'
                      + "scan_num"
                      + '\t'
                      + "quantification_peak_sum"
                      + '\t'
                      + "quantification_peak_area"
                      + '\t'
                      + "svr_score"
                      + '\t'
                      + "quality_score"
                      + '\t'
                      + "mzs,rts,ints"
                      + '\n';
      printWriter.print(header);

      for (int i = 0; i < scorePepAll.size(); i++) {
        int size = scorePepAll.get(i).size();
        List<XIC> featureIsotopes = trailsAll.get(i);
        double peak_sum = 0, peak_area = 0;
        for (XIC isotope: featureIsotopes) {
          peak_sum += isotope.getPeakSum();
          peak_area += isotope.getPeakArea();
        }

        for (int j = 0; j < size; j++) {
          int id = i + 1;
          XIC trail = featureIsotopes.get(j);
          int scannum = trail.getScanNum();
          ArrayList<Double> ints = trail.getIntensities();
          ArrayList<Double> mzs = trail.getMassChargeRatios();
          ArrayList<Double> rts = trail.getRetentionTimes();

          String mzsStr =  mzs.get(0).toString();
          String rtsStr =  rts.get(0).toString();
          String intsStr =  ints.get(0).toString();
          for (int k = 1; k < scannum; k++) {
            mzsStr += '\t' + mzs.get(k).toString();
            rtsStr += '\t' + rts.get(k).toString();
            intsStr += '\t' + ints.get(k).toString();
          }

          String data = "";
          data += id;
          data += '\t' + scorePepAll.get(i).get(j).getR().toString();
          data += '\t' + scorePepAll.get(i).get(j).getL().toString();
          data += '\t' + scoreZVal.get(i).toString();
          data += '\t' + scoreIsoNum.get(i).toString();
          data += '\t' + scoreIntShape.get(i).toString();
          data += '\t' + scoreIsoDistr.get(i).toString();
          data += '\t' + intensityPercentage.get(i).toString();
          data += '\t' + Integer.toString(scannum);
          data += '\t' + Double.toString(peak_sum);
          data += '\t' + Double.toString(peak_area);
          data += '\t' + "0"; // placeholder for svr_score
          data += '\t' + "0"; // placeholder for quality score
          data += '\t' + mzsStr;
          data += '\t' + rtsStr;
          data += '\t' + intsStr;
          data += '\n';
          printWriter.print(data);
        }
      }
      // closing writer connection
      printWriter.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }


  /**
   * Search for isotopes of all charge states based on raw signal data before local maximal intensities;
   * calculate
   *
   * @param window_size the max distance to search isotopes
   * @param num_of_z z values: 1 ... z
   */
  //@author JCZhong
  public void searchIsotopeWithSignalToMS1(MSOneTrailSet spec, int window_size, int num_of_z) {
    //Check all raw mz value, find all possible isotope signal in raw data
    //for 1 charge, mz + 1 with ppm
    //for 2 charge, mz + 0.5 with ppm
    //for N charge, mz + 1/N with ppm
    // Search max near one m/z vertically[current mass value,+N(N=10)] default charge value is 10
    int id=1;


//    List<SignalIsotope> liSigIsotope = new ArrayList<>();
    int num_of_Isotope = 8;//
    for (int width = 0; width < IDX; width++) {
      int scanMassCount = mzMap[width].length;
      for (int height = 0; height < scanMassCount-1; height++) {
        SignalIsotope sIsotope = null;
        int iCurPos = height;
        for (int iZ = num_of_z; iZ >= 1; iZ--) {
          List<Integer> listZIsotopePos = null;
          List<Double> listZIsotopeValue = null;

          int iIsotpe = 1;
          boolean bzHighChargeCover = false;
          boolean bzFindNew = false;

          double dHighMzofIsotopeWithZ = mzMap[width][height] * (1 + ppm) + num_of_Isotope * (C13 - C12) / iZ;

          int iNextHighIsotopePos = getNearpos(mzMap[width], dHighMzofIsotopeWithZ, iCurPos, scanMassCount);//the maximum pos
          if (iNextHighIsotopePos + 1 <= scanMassCount) {
            iNextHighIsotopePos++;//the pos  greater than the dHighMzofIsotopeWithZ
          }

          int iCurrentNearPos = 0;
          int iCurrentZ_num_Isotope =0;
          for (; iIsotpe <= num_of_Isotope; iIsotpe++) {//find the nearest pos
            iCurrentNearPos = 0;
            bzHighChargeCover = false;
            //IF Z is half of num_of_z,we should j
            if (iZ <= num_of_z / 2) {
              int tmpZ = iZ;
              int tmpIsotope = iIsotpe;
              while (tmpZ+iZ <= num_of_z) {

                tmpZ += iZ;
                tmpIsotope = iIsotpe * tmpZ / iZ;

                if (sIsotope!= null && sIsotope.lliIstopePos!=null && sIsotope.lliIstopePos.get(tmpZ) != null) {
                  if (sIsotope.lliIstopePos.get(tmpZ).size() >= tmpIsotope) {
                    bzHighChargeCover = true;
                    if (listZIsotopePos == null) {
                      listZIsotopePos = new ArrayList<>();
                      listZIsotopeValue = new ArrayList<>();
                    }
                    listZIsotopePos.add(sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1));
                    listZIsotopeValue.add(mzMap[width][sIsotope.lliIstopePos.get(tmpZ).get(tmpIsotope-1)]);
                  }
                }
              }
            }
            if (!bzHighChargeCover) {//没有覆盖
              double dRealMzOfIsotopeWithZ = mzMap[width][height];
              int iNextIsotopePos = getNearpos(mzMap[width], dRealMzOfIsotopeWithZ + iIsotpe * (C13 - C12) / iZ, iCurPos + 1, iNextHighIsotopePos);
              iCurrentNearPos = isInPPM(mzMap[width], dRealMzOfIsotopeWithZ, iNextIsotopePos, iIsotpe, iZ);
            }
            if (iCurrentNearPos > 0) {

              if (sIsotope == null) {
                sIsotope = new SignalIsotope();
                sIsotope.iRTpos = width;
                sIsotope.iMzPos = height;
                sIsotope.mzvalue = mzMap[width][height];
                sIsotope.lliIstopePos = new HashMap<>();
                sIsotope.lliIstopeValue = new HashMap<>();
              }
              if (listZIsotopePos == null) {
                listZIsotopePos = new ArrayList<>();
                listZIsotopeValue = new ArrayList<>();
              }
              listZIsotopePos.add(iCurrentNearPos);
              listZIsotopeValue.add(mzMap[width][iCurrentNearPos]);
              iCurrentZ_num_Isotope = iIsotpe;
              bzFindNew = true;

            } else if(!bzHighChargeCover) {
              break;
            }

          }
          if (bzFindNew)//NO COVERED BY HIGH CHARGE AND FIND NEW
          {
            sIsotope.lliIstopePos.put(iZ,listZIsotopePos);
            sIsotope.lliIstopeValue.put(iZ,listZIsotopeValue);
            spec.readFromeTrail(mzMap[width][height],RT[width],iZ,
                    iCurrentZ_num_Isotope,width, intensityMap[width][height],true);

          }

        }


      }
    }
  }
}