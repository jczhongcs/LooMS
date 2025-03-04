package edu.uw.waterlooms.ms1;

import edu.uw.waterlooms.entity.Pair;
import edu.uw.waterlooms.entity.SVRScore;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.peptideMatch.MSOneTrailSet;
import edu.uw.waterlooms.service.ParameterService;

import java.io.*;
import java.util.*;

public class FeatureSelect {
  private int id_index;
  private int mz_index;
  private int rt_index;
  private int z_index;
  private int isonum_index;
  private int int_shape_index;
  private int iso_distr_index;
  private int intensity_area_percentage_index;
  private int scan_num_index;
  private int quantification_peaks_sum_index;
  private int quantification_peaks_area_index;
  private int svr_index;
  private int quality_index;
  private int invalidVal;
  private double mz_error;
  private double rt_error;
  private List<Double[]> model_result; // read _svr_nn_score
  private ArrayList<XIC> precursor_trails;

  public MSOneTrailSet getSpec() {
    return spec;
  }


  private MSOneTrailSet spec;

  enum ML {
    SVR,
    NN
  }
  /**
   * Select charge state based on SVR score; one charge state of a feature
   *
   * @param oldFilePath file path with defined suffix
   * @throws IOException
   */
  public void selectFeature(String oldFilePath) throws IOException {
    setParameters();
    ML ml = ML.SVR;
    String filepath = oldFilePath.replaceFirst("[.][^.]+$", "");
    model_result = readFile(filepath + "_SVRScore.tsv");
    selectCharge();
//    model_result = removeIso(model_result);
    writeFile(filepath + "_featureOneZ.tsv", model_result, ml);
  }


  /**
   * Generate the final feature list; last step
   *
   * @param oldFilePath file path without defined suffix TODO: Refactor this
   * @param workingDirectory String containing the working directory
   * @param rawFileName String containing the name of the MSConvert mzXML file (converted from RAW)
   * @throws IOException
   */
  public void finalizeFeatureToMS1(String oldFilePath, String workingDirectory, String rawFileName) throws IOException {
    setParameters();
    ML ml = ML.NN;
    String filepath = oldFilePath.replaceFirst("[.][^.]+$", "");
    model_result = readFile(filepath + "_NNScore.tsv");
    model_result = clusterFeature(model_result, mz_error, rt_error);

//    List<SVRScore> svrScores = new ArrayList<>();
//    for (Double[] element : model_result) {
//      SVRScore ms1Precursor =
//              new SVRScore(
//                      element[0],
//                      element[1],
//                      element[2],
//                      element[3],
//                      element[4],
//                      element[5],
//                      element[6],
//                      element[7],
//                      element[8],
//                      element[9],
//                      element[10],
//                      element[11],
//                      element[12],
//                      element[13],
//                      element[14]);
//      svrScores.add(ms1Precursor);
//    }
//
//    // TODO: For each model_result row, save them into svrScores.
//
//    // Take trails from precursor isotopes
//    abstractPrecursorTrails();
//    writeFile(filepath + "_prefinalized_precursors.tsv", model_result, ml);

    // Finalize precursor file
    model_result = removeIso(model_result);
    model_result.sort(Comparator.comparing(l -> l[quality_index]));
//    Collections.reverse(model_result);
//    writeFinalFile(filepath + "_precursors.tsv", model_result);
    toMS1WithFeature(model_result);
//    return svrScores;
  }



  /**
   * Generate the final feature list; last step
   *
   * @param oldFilePath file path without defined suffix TODO: Refactor this
   * @param workingDirectory String containing the working directory
   * @param rawFileName String containing the name of the MSConvert mzXML file (converted from RAW)
   * @throws IOException
   */
  public List<SVRScore> finalizeFeature(String oldFilePath, String workingDirectory, String rawFileName) throws IOException {
    setParameters();
    ML ml = ML.NN;
    String filepath = oldFilePath.replaceFirst("[.][^.]+$", "");
    model_result = readFile(filepath + "_NNScore.tsv");
    model_result = clusterFeature(model_result, mz_error, rt_error);

    List<SVRScore> svrScores = new ArrayList<>();
    for (Double[] element : model_result) {
      SVRScore ms1Precursor =
          new SVRScore(
              element[0],
              element[1],
              element[2],
              element[3],
              element[4],
              element[5],
              element[6],
              element[7],
              element[8],
              element[9],
              element[10],
              element[11],
              element[12],
              element[13],
              element[14]);
      svrScores.add(ms1Precursor);
    }

    // TODO: For each model_result row, save them into svrScores.

    // Take trails from precursor isotopes
    abstractPrecursorTrails();
    writeFile(filepath + "_prefinalized_precursors.tsv", model_result, ml);

    // Finalize precursor file
    model_result = removeIso(model_result);
    model_result.sort(Comparator.comparing(l -> l[quality_index]));
    Collections.reverse(model_result);
    writeFinalFile(filepath + "_precursors.tsv", model_result);
    return svrScores;
  }

  /**
   * Helper function for constructor. Sets parameter values defined above.
   *
   */
  private void setParameters() {
    id_index = ParameterService.getIdIndex();
    mz_index = ParameterService.getMzIndex();
    rt_index = ParameterService.getRtIndex();
    z_index = ParameterService.getzIndex();
    isonum_index = ParameterService.getIsonumIndex();
    int_shape_index = ParameterService.getIntShapeIndex();
    iso_distr_index = ParameterService.getIsoDistrIndex();
    intensity_area_percentage_index = ParameterService.getIntensityAreaPercentageIndex();
    scan_num_index = ParameterService.getScanSumIndex();
    quantification_peaks_sum_index = ParameterService.getPeaksSumIndex();
    quantification_peaks_area_index = ParameterService.getPeaksAreaIndex();
    svr_index = ParameterService.getSvrIndex();
    quality_index = ParameterService.getQualityIndex();
    invalidVal = ParameterService.getInvalidVal();
    mz_error = ParameterService.getMzError();
    rt_error = ParameterService.getRtError();
  }

  public void selectCharge() {
    model_result.sort(Comparator.comparing(l -> l[rt_index]));
    model_result.sort(Comparator.comparing(l -> l[mz_index]));
    for (int i = 0; i < model_result.size(); ) {
      double cur_mz = model_result.get(i)[mz_index];
      double cur_rt = model_result.get(i)[rt_index];
      if (i == model_result.size() - 1) {
        break;
      }
      int best_index = i;

      double best_score = model_result.get(i)[svr_index];
      for (int j = 1; i + j < model_result.size(); j++) {
        double next_mz = model_result.get(i + j)[mz_index];
        double next_rt = model_result.get(i + j)[rt_index];
        double next_svr = model_result.get(i + j)[svr_index];
        if (cur_mz == next_mz && cur_rt == next_rt) {
          if (next_svr > best_score) {
            model_result.get(best_index)[svr_index] = (double) invalidVal;
            best_score = next_svr;
            best_index = i + j;
          } else {
            model_result.get(i + j)[svr_index] = (double) invalidVal;
          }
        } else {
          i += j;
          break;
        }
        if (i + j == model_result.size() - 1) {
          i = model_result.size();
        }
      }
    }
    model_result.sort(Comparator.comparing(l -> l[id_index]));
    int total = 0;
    for (int i = 0; i < model_result.size(); ) {
      int curPepNum = model_result.get(i)[isonum_index].intValue();
      int fst = 0;
      int last = curPepNum - 1;
      boolean is_valid = true;
      while (fst < curPepNum && model_result.get(i + fst)[svr_index] == invalidVal) {
        fst++;
      }
      while (last >= 0 && model_result.get(i + last)[svr_index] == invalidVal) {
        last--;
      }
      if (fst >= last) { // if there is only 1 left in the group or all are not left
        is_valid = false;
      } else {
        for (int j = fst; j <= last; j++) {
          if (model_result.get(i + j)[svr_index] == invalidVal) {
            is_valid = false;
            break;
          }
        }
      }
      if (!is_valid) {
        for (int j = 0; j < curPepNum; j++) {
          model_result.get(i + j)[svr_index] = (double) invalidVal;
        }
      } else {
        total++;
      }
      i += curPepNum;
    }
    // Delete invalid score
    List<Double[]> valid_model_result = new ArrayList<>();
    for (int i = 0; i < model_result.size(); i++) {
      if (model_result.get(i)[svr_index] != (double) invalidVal) {
        valid_model_result.add(model_result.get(i));
      }
    }
    model_result = valid_model_result;
    // Change isotope number in a group
    for (int i = 0; i < model_result.size(); ) {
      int cur_pepnum = 1;
      int cur_id = model_result.get(i)[id_index].intValue();
      int j = i + 1;
      for (; j < model_result.size(); j++) {
        int next_id = model_result.get(j)[id_index].intValue();
        if (cur_id == next_id) {
          cur_pepnum++;
        } else {
          break;
        }
      }
      for (int k = 0; k < cur_pepnum; k++) {
        model_result.get(i + k)[isonum_index] = (double) cur_pepnum;
      }
      i = j;
    }
    System.out.println("selectCharge completed");
  }
  /**
   * Cluster similar features
   *
   * @param isoList input isotopic feature lists
   * @param mz_error mz error for clustering
   * @param rt_error rt errot for clustering
   * @return clustered feature list
   */
  public List<Double[]> clusterFeature(List<Double[]> isoList, double mz_error, double rt_error) {
    List<Double[]> featureGroup = new ArrayList<>();
    List<Double[]> isoAllFeature = removeIso(isoList);
    isoAllFeature.sort(Comparator.comparing(l -> l[rt_index]));
    isoAllFeature.sort(Comparator.comparing(l -> l[z_index]));
    isoAllFeature.sort(Comparator.comparing(l -> l[mz_index]));
    for (int i = 0; i < isoAllFeature.size();) {
      double mz_iso1 = isoAllFeature.get(i)[mz_index];
      double rt_iso1 = isoAllFeature.get(i)[rt_index];
      double z_iso1 = isoAllFeature.get(i)[z_index];
      int isonum_iso1 = isoAllFeature.get(i)[isonum_index].intValue();
      double cur_hi_score = isoAllFeature.get(i)[quality_index];

      List<Double[]> newIso = new ArrayList<>();
      newIso.add(isoAllFeature.get(i).clone());
      int j = i + isonum_iso1;
      for (; j < isoAllFeature.size(); ) {
        double mz_iso2 = isoAllFeature.get(j)[mz_index];
        double rt_iso2 = isoAllFeature.get(j)[rt_index];
        double z_iso2 = isoAllFeature.get(j)[z_index];
        double score_iso2 = isoAllFeature.get(j)[quality_index];
        if (mz_iso2 <= mz_iso1 * (1 + mz_error)
                && mz_iso2 >= mz_iso1 * (1 - mz_error)
                && rt_iso2 >= rt_iso1 - rt_error
                && rt_iso2 <= rt_iso1 + rt_error
                && z_iso1 == z_iso2) {
          if (score_iso2 > cur_hi_score) {
            cur_hi_score = score_iso2;
            newIso.add(isoAllFeature.get(j).clone());
          }
          isoAllFeature.remove(j);
        } else if (mz_iso2 >= mz_iso1 * (1 + mz_error)) {
          break;
        } else if (mz_iso2 <= mz_iso1 * (1 - mz_error)) {
          break;
        } else {
          j++;
        }
      }
      featureGroup.addAll(newIso);
      i++;
    }

    // Add isotopes back
    featureGroup.sort(Comparator.comparing(l -> l[id_index]));
    int j = 0;
    for (int i = 0; i < featureGroup.size();) {
      Double[] feature = featureGroup.get(i);
      for (; j < isoList.size();) {
        int isonum = isoList.get(j)[isonum_index].intValue();
        if (feature[id_index].intValue() == isoList.get(j)[id_index].intValue()) {
          for (int k = 1; k < isonum; k++) {
            featureGroup.add(i + k, isoList.get(j + k));
          }
          i += isonum;
          j += isonum;
          break;
        } else {
          j += isonum;
        }
      }
    }
    System.out.println("Peptide # after deleting duplicates: " + featureGroup.size());
    System.out.println("clusterFeature completed");
    return featureGroup;
  }

  /**
   * Removing redundant information of group
   *
   * @param isolist input isotopic feature lists
   * @return list leaving with one main isotope of the group
   */
  public List<Double[]> removeIso(List<Double[]> isolist) {
    List<Double[]> list = new ArrayList<>();
    for (int i = 0; i < isolist.size(); ) {
      list.add(isolist.get(i));
      i += isolist.get(i)[isonum_index];
    }
    return list;
  }

  public void abstractPrecursorTrails() {
    List<Double[]> list = model_result;
    precursor_trails = new ArrayList<>();

    for (int i = 0; i < list.size(); i++) {
      int scannum = list.get(i)[scan_num_index].intValue();
      ArrayList<Double> ints = new ArrayList<>(), mzs = new ArrayList<>(), rts = new ArrayList<>();
      for (int k = 0; k < scannum; k++) {
        mzs.add(list.get(i)[quality_index + 1 + k]);
        rts.add(list.get(i)[quality_index + 1 + scannum + k]);
        ints.add(list.get(i)[quality_index + 1 + scannum * 2 + k]);
      }
      XIC trail = new XIC(ints, mzs, rts, list.get(i)[quantification_peaks_sum_index], list.get(i)[quantification_peaks_area_index]);
      precursor_trails.add(trail);
    }
  }

  public ArrayList<XIC> getPrecursorTrails() {
    return precursor_trails;
  }

  static List<Double[]> readFile(String filename) throws IOException {
    List<Double[]> result = new ArrayList<>();
    FileReader fileReader = new FileReader(filename);
    BufferedReader bufferedReader = new BufferedReader(fileReader);
    String line;
    boolean IsfirstLine = true;
    while ((line = bufferedReader.readLine()) != null) {
      // process first line
      if (IsfirstLine == true) {
        IsfirstLine = false;
      } else {
        // process line
        String[] words = line.split("\\t");
        Double[] nums = new Double[words.length];
        for (int i = 0; i < words.length; i++) {
          nums[i] = Double.valueOf(words[i]);
        }
        result.add(nums);
      }
    }
    return result;
  }

  void writeFile(String filename, List<Double[]> list, ML ml) throws IOException {
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
                      + "quantification_peaks_sum"
                      + '\t'
                      + "quantification_peaks_area"
                      + '\t'
                      + "svr_score"
                      + '\t'
                      + "quality_score"
                      + '\t'
                      + "mzs,rts,ints"
                      + '\n';
      printWriter.print(header);
      for (int i = 0; i < list.size(); i++) {
        Integer z = list.get(i)[z_index].intValue();
        Integer iso_num = list.get(i)[isonum_index].intValue();
        int scannum = list.get(i)[scan_num_index].intValue();
        String mzsStr = list.get(i)[quality_index + 1].toString();
        String rtsStr =  list.get(i)[quality_index + 1 + scannum].toString();
        String intsStr = list.get(i)[quality_index + 1 + scannum * 2].toString();
        for (int k = 1; k < scannum; k++) {
          mzsStr += "\t" + list.get(i)[quality_index + 1 + k];
          rtsStr += "\t" + list.get(i)[quality_index + 1 + scannum + k];
          intsStr += "\t" + list.get(i)[quality_index + 1 + scannum * 2 + k];
        }

        String data = "";
        data += list.get(i)[id_index].toString();
        data += '\t' + list.get(i)[mz_index].toString();
        data += '\t' + list.get(i)[rt_index].toString();
        data += '\t' + z.toString();
        data += '\t' + iso_num.toString();
        data += '\t' + list.get(i)[int_shape_index].toString();
        data += '\t' + list.get(i)[iso_distr_index].toString();
        data += '\t' + list.get(i)[intensity_area_percentage_index].toString();
        data += '\t' + Integer.toString(scannum);
        data += '\t' + list.get(i)[quantification_peaks_sum_index].toString();
        data += '\t' + list.get(i)[quantification_peaks_area_index].toString();
        data += '\t' + list.get(i)[svr_index].toString();
        if (ml == ML.NN) {
          data += '\t' + list.get(i)[quality_index].toString();
        } else if (ml == ML.SVR) {
          data += '\t' + "0";
        }
        data += '\t' + mzsStr;
        data += '\t' + rtsStr;
        data += '\t' + intsStr;
        data += '\n';
        printWriter.print(data);
      }
      // closing writer connection
      printWriter.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  void writeFinalFile(String filename, List<Double[]> list) throws IOException {
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
                      + "quantification_peaks_sum"
                      + '\t'
                      + "quantification_peaks_area"
                      + '\t'
                      + "svr_score"
                      + '\t'
                      + "quality_score"
                      + '\t'
                      + "mzs"
                      + '\t'
                      + "rts"
                      + '\t'
                      + "ints"
                      + '\n';
      printWriter.print(header);
      Integer id = 1;
      for (int i = 0; i < list.size(); i++) {
        Integer z = list.get(i)[z_index].intValue();
        Integer iso_num = list.get(i)[isonum_index].intValue();
        int scannum = list.get(i)[scan_num_index].intValue();
        String mzsStr = "[" + list.get(i)[quality_index + 1].toString();
        String rtsStr = "[" + list.get(i)[quality_index + 1 + scannum].toString();
        String intsStr = "[" + list.get(i)[quality_index + 1 + scannum * 2].toString();
        for (int k = 1; k < scannum; k++) {
          mzsStr += "," + list.get(i)[quality_index + 1 + k];
          rtsStr += "," + list.get(i)[quality_index + 1 + scannum + k];
          intsStr += "," + list.get(i)[quality_index + 1 + scannum * 2 + k];
        }
        mzsStr += "]";
        rtsStr += "]";
        intsStr += "]";

        String data = "";
        data += id.toString();
        data += '\t' + list.get(i)[mz_index].toString();
        data += '\t' + list.get(i)[rt_index].toString();
        data += '\t' + z.toString();
        data += '\t' + iso_num.toString();
        data += '\t' + list.get(i)[int_shape_index].toString();
        data += '\t' + list.get(i)[iso_distr_index].toString();
        data += '\t' + list.get(i)[intensity_area_percentage_index].toString();
        data += '\t' + Integer.toString(scannum);
        data += '\t' + list.get(i)[quantification_peaks_sum_index].toString();
        data += '\t' + list.get(i)[quantification_peaks_area_index].toString();
        data += '\t' + list.get(i)[svr_index].toString();
        data += '\t' + list.get(i)[quality_index].toString();
        data += '\t' + mzsStr;
        data += '\t' + rtsStr;
        data += '\t' + intsStr;
        data += '\n';
        printWriter.print(data);
        id++;
      }
      // closing writer connection
      printWriter.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  void toMS1WithFeature( List<Double[]> list) throws IOException {
    if(spec==null) {
      spec = new MSOneTrailSet();
    }
    spec.readFromeListTrail(list,false);
    spec.generateRTSandMapRT();

  }
}
