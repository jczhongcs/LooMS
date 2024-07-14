package edu.uw.waterlooms.ms2;

import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.msutil.ScanEntry;
import java.io.IOException;

import java.util.*;
import java.util.stream.Collectors;

class c_min_mz implements Comparator<ScanEntry> {
  // Used for sorting in ascending order of
  // roll number
  public int compare(ScanEntry a, ScanEntry b)
  {
    if (a.mz < b.mz) {
      return -1;
    } else {
      return 1;
    }
  }
}


// *** this class holds extracted MS2 feature Extraction ***
// *** final output:
//        extracted XICs by peak points from 1 isolation window
// *** functions including:
//        Finding local max points (peaks) in scans from 1 isolation window
//        filtering local max points
//        feature extraction
public class MS2Extraction {
  ArrayList<ArrayList<ScanEntry>> isoWin;
  private ArrayList<ScanEntry> local_maxs = new ArrayList<>();
  private ArrayList<XIC> xics = new ArrayList<>();
  Map<ScanEntry, XIC> MS2_xics = new HashMap<>();
  // Local max criteria

  // TODO: Pass in MS2 XIC Parameter File
  // TODO: Could just specify this in FeatureDetect file
  private static final int MS2TorPPM = 7; // TODO: check this value
  private static final int MS2TorScans = 1;
  private static final double PPMCalNumU = 1 + MS2TorPPM  * Math.pow(10, -6);
  private static final double PPMCalNumL = 1 - MS2TorPPM  * Math.pow(10, -6);
  // feature extraction criteria
//    private static double MS2IntCutoff = 0.1; // TODO: check this value
  private static final double MS2TorRt = 0.5;

  public MS2Extraction(ArrayList<ArrayList<ScanEntry>> isoWin) throws IOException {
    this.isoWin = isoWin;
    find_local_max(isoWin);
    feature_extraction();
    filter_local_max();
  }

  public Map<ScanEntry, XIC> getMS2_xics() {
    return MS2_xics;
  }

  private void find_local_max(ArrayList<ArrayList<ScanEntry>> isoWin) throws IOException {
    int num_specs = isoWin.size();
    for (int j = 0; j < num_specs; ++ j) {
      ArrayList<ScanEntry> spectrum = isoWin.get(j);
      for (ScanEntry point : spectrum) {
        double mz_base = point.mz;
        double int_base = point.intensity;
        double rt_base = point.rt;
        double MZL = mz_base * PPMCalNumL;
        double MZU = mz_base * PPMCalNumU;
        int j_start = Math.max(j - MS2TorScans, 0);
        int j_end = Math.min(j + MS2TorScans, num_specs - 1);
        if (j_end - j_start >= 2) {
          boolean has_trace = true;
          for (int j_c = j_start; j_c <= j_end; ++j_c) {
            ArrayList<ScanEntry> spec_comp = isoWin.get(j_c);
            int index = Collections.binarySearch(spec_comp, new ScanEntry(MZL, 0, 0), new c_min_mz());
            if (index < 0) {
              index = -(index + 1);
            }
            int indexU = Collections.binarySearch(spec_comp, new ScanEntry(MZU, 0, 0), new c_min_mz());
            if (indexU < 0) {
              indexU = -(indexU + 1);
            }
            if (index == indexU) {
              has_trace = false;
              break;
            }
            for (int q = index; q < indexU; ++q) {
              double cur_mz = spec_comp.get(q).mz;
              double cur_int = spec_comp.get(q).intensity;
              if (cur_mz > MZU) break;
              if (cur_int > int_base) {
                has_trace = false;
                break;
              }
            }
          }
          if (has_trace) {
            local_maxs.add(new ScanEntry(mz_base, rt_base, int_base));
          }
        }
      }
    }

    /* TODO: Refactor this to accept a parameter
    String out_path = "";
    FileWriter fstream_t = new FileWriter(out_path + "local_max.mgf");
    BufferedWriter info_t = new BufferedWriter(fstream_t);
    Collections.sort(local_maxs);
    for (ScanEntry lm: local_maxs) {
      info_t.write(lm.mz + " " + lm.rt);
      info_t.newLine();
    }
    info_t.close(); */
  }

  private void filter_local_max() {
    int index = 0;
    for (XIC xic_cur: xics) {
      List<Double> ints = xic_cur.getIntensities();
      List<Double> rts = xic_cur.getRetentionTimes();
      List<Double> mzs = xic_cur.getMassChargeRatios();
      // find local max in extracted xics
      List<ScanEntry> local_maxs_e = new ArrayList<>();
      int num_points = ints.size();
      for (int q = MS2TorScans; q < num_points - MS2TorScans; ++q) {
        int start_index = q - MS2TorScans;
        int end_index = q + MS2TorScans;
        double cur_int = ints.get(q);
        double cur_mz = mzs.get(q);
        double cur_rt = rts.get(q);
        boolean is_max = true;
        for (int r = start_index; r <= end_index; ++r) {
          double next_int = ints.get(r);
          if (cur_int < next_int) {
            is_max = false;
            break;
          }
        }
        if (is_max) {
          local_maxs_e.add(new ScanEntry(cur_mz, cur_rt, cur_int));
        }
      }
      // trim xic
      int right_bound = ints.size();
      int left_bound = 0;
      // trim between local maxs
      for (int i = 0; i < local_maxs_e.size() - 1; ++ i) {
        ScanEntry cur_local_max = local_maxs_e.get(i);
        ScanEntry next_local_max = local_maxs_e.get(i + 1);
        int cur_index = ints.indexOf(cur_local_max.intensity);
        int next_index = ints.indexOf(next_local_max.intensity);
        boolean indep = false;
        for (int q = cur_index + 1; q < next_index; ++q) {
          double cur_int = ints.get(q);
          if (cur_int <= cur_local_max.intensity * 0.5) {
            right_bound = q + 1;
            indep = true;
            break;
          }
        }
        if (indep) {
          add_to_map(ints, rts, mzs, left_bound, right_bound, cur_local_max.rt, cur_local_max.intensity);
          left_bound = right_bound;
          right_bound = ints.size();
        }
      }
      ArrayList<Double> ints_f = new ArrayList<>(ints.subList(left_bound, right_bound));
      ArrayList<Double> rts_f = new ArrayList<>(rts.subList(left_bound, right_bound));
      int max_index = ints_f.indexOf(Collections.max(ints_f));
      add_to_map(ints, rts, mzs, left_bound, right_bound, rts_f.get(max_index), ints_f.get(max_index));
      ++ index;
    }
  }

  private void feature_extraction() {
    // Extract MS2 XIC and refine for each detected local maximum
    for (ScanEntry local_max: local_maxs) {
      double rt = local_max.rt;
      double mz = local_max.mz;
      double MZL = mz * PPMCalNumL;
      double MZU = mz * PPMCalNumU;
      // Extract XIC
      List<ArrayList<ScanEntry>> scan_f = isoWin.stream()
          .filter(s -> s.get(0).rt <= rt + MS2TorRt && s.get(0).rt >= rt - MS2TorRt)
          .collect(Collectors.toList());
      int num_scans = scan_f.size();
      ArrayList<Double> rts_left = new ArrayList<>();
      ArrayList<Double> mzs_left = new ArrayList<>();
      ArrayList<Double> ints_left = new ArrayList<>();
      int rt_index = scan_f.stream().map(s -> s.get(0).rt).collect(Collectors.toList()).indexOf(rt);
      for (int m = rt_index; m >= 0; -- m) {
        ArrayList<ScanEntry> spec_compare = scan_f.get(m);
        int index = Collections.binarySearch(spec_compare, new ScanEntry(MZL, 0, 0), new c_min_mz());
        if (index < 0) {
          index = -(index + 1);
        }
        if (index >= spec_compare.size()) {
          break;
        }
        double cur_mz = spec_compare.get(index).mz;
        double cur_int = spec_compare.get(index).intensity;
        double cur_rt = spec_compare.get(index).rt;
        if (cur_mz > MZU) break;
        ints_left.add(0, cur_int);
        mzs_left.add(0, cur_mz);
        rts_left.add(0, cur_rt);
      }
      for (int m = rt_index + 1; m < num_scans; ++ m) {
        ArrayList<ScanEntry> spec_compare = scan_f.get(m);
        int index = Collections.binarySearch(spec_compare, new ScanEntry(MZL, 0, 0), new c_min_mz());
        if (index < 0) {
          index = -(index + 1);
        }
        if (index >= spec_compare.size()) {
          break;
        }
        double cur_mz = spec_compare.get(index).mz;
        double cur_int = spec_compare.get(index).intensity;
        double cur_rt = spec_compare.get(index).rt;
        if (cur_mz > MZU) break;
        ints_left.add(cur_int);
        mzs_left.add(cur_mz);
        rts_left.add(cur_rt);
      }
      XIC xic = new XIC(ints_left, mzs_left, rts_left);
      xics.add(xic);
    }
  }

  private static double dot_product(List<Double> x, List<Double> y, int len) {
    double dot_prod = 0;
    for (int d = 0; d < len; ++ d) {
      dot_prod += x.get(d) * y.get(d);
    }
    return dot_prod;
  }

  private void add_to_map(List<Double> ints, List<Double> rts, List<Double> mzs, int left_bound, int right_bound, double rt, double intensity) {
    // Sanity check to terminate execution if the trimming somehow got the boundaries wrong
    if (left_bound >= right_bound){
      return;
    }
    //System.out.println("left_bound : right_bound: [" + left_bound + ":" + right_bound + "]");
    ArrayList<Double> ints_f = new ArrayList<>(ints.subList(left_bound, right_bound));
    ArrayList<Double> rts_f = new ArrayList<>(rts.subList(left_bound, right_bound));
    ArrayList<Double> mzs_f = new ArrayList<>(mzs.subList(left_bound, right_bound));
    double centroid_mz = dot_product(ints_f, mzs_f, ints_f.size()) / ints_f.stream().mapToDouble(a -> a).sum();
    centroid_mz = Math.round(centroid_mz * 100000d) / 100000d;
    // add current XIC to dictionary
    XIC xic = new XIC(ints_f, mzs_f, rts_f);
    ScanEntry local_max = new ScanEntry(centroid_mz, rt, intensity);
    if (!MS2_xics.containsKey(local_max)) {
      MS2_xics.put(local_max, xic);
    }
  }
}
