package edu.uw.waterlooms.service;

// TODO: This class is depreciated until the refactor occurs with Veronica's work.

import edu.uw.waterlooms.entity.Precursor;
import edu.uw.waterlooms.msutil.ScanEntry;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.msutil.SpecEntry;

import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

// *** This class holds matching algorithms for precursor and product ion ***
public class PrecursorMatch {
  private Map<ScanEntry, XIC> MS2_xics;
  private ArrayList<Map.Entry<Double, Precursor>> MS1_xics;
  Map<Precursor, ArrayList<SpecEntry>> spectra = new HashMap<>();
  Map<Precursor, ArrayList<SpecEntry>> MS1_2_cor = new HashMap<>();
  Map<SpecEntry, ArrayList<Precursor>> MS2_1_cor = new HashMap<>();
  BufferedWriter info;
  // matching criteria
  private static double XICMatchTorRt = 0.6;
  private static double XICMatchScore = 0.9;
  private static double RtInter = 1.0 / 150;

  public PrecursorMatch(Map<ScanEntry, XIC> MS2_xics, ArrayList<Map.Entry<Double, Precursor>> MS1_xics, BufferedWriter info) throws IOException {
    this.MS2_xics = MS2_xics;
    this.MS1_xics = MS1_xics;
    this.info = info;
    get_matching_score();
    generate_spec();
  }

  private void get_matching_score() {
    // Double is score for PQ, XIC is for MS2
    ArrayList<Precursor> ps = new ArrayList<>();
    for (Map.Entry<ScanEntry, XIC> entry : MS2_xics.entrySet()) {
      ScanEntry local_max = entry.getKey();
      XIC MS2_xic = entry.getValue();
      double rt = local_max.rt;
      double mz = local_max.mz;
      List<Double> MS2_xic_rt = MS2_xic.getRetentionTimes();
      List<Double> MS2_xic_ints = MS2_xic.getIntensities();
      List<Double> MS2_xic_rt_i = new ArrayList<>();
      List<Double> MS2_xic_ints_i = new ArrayList<>();
      if (MS2_xic_ints.size() <= 2) continue;
      // interpolate to same interval size
      interpolate(MS2_xic_rt, MS2_xic_ints, MS2_xic_ints_i, MS2_xic_rt_i);
      int ms2_len = MS2_xic_ints_i.size();
      // get all precursor ions in the rt range
      List<Map.Entry<Double, Precursor>> MS1_xic_f = MS1_xics.stream()
          .filter(x -> x.getKey() >= rt - XICMatchTorRt && x.getKey() <= rt + XICMatchTorRt)
          .collect(Collectors.toList());

      // get matching scores for each MS1 xic in range and add to PQ
      for (Map.Entry<Double, Precursor> xic_f: MS1_xic_f) {
        Precursor MS1_xic = xic_f.getValue();
        List<Double> MS1_xic_rt = MS1_xic.xic_rt;
        List<Double> MS1_xic_ints = MS1_xic.xic_ints;
        List<Double> MS1_xic_rt_i = new ArrayList<>();
        List<Double> MS1_xic_ints_i = new ArrayList<>();
        if (MS1_xic_rt.size() <= 2) {
          if (!ps.contains(MS1_xic)) {
            ps.add(MS1_xic);
          }
          continue;
        }
        // linear interpolate to same interval size
        interpolate(MS1_xic_rt, MS1_xic_ints, MS1_xic_ints_i, MS1_xic_rt_i);
        int ms1_len = MS1_xic_ints_i.size();
        int max_rti_1 = MS1_xic_ints_i.indexOf(Collections.max(MS1_xic_ints_i));
        int max_rti_2 = MS2_xic_ints_i.indexOf(Collections.max(MS2_xic_ints_i));
        // convert to same length
        int left_size = Math.min(max_rti_1, max_rti_2);
        int right_size = Math.min(ms1_len - max_rti_1, ms2_len - max_rti_2) - 1;
        List<Double> MS1_xic_rt_t = new ArrayList<>();
        List<Double> MS1_xic_ints_t = new ArrayList<>();
        trim(MS1_xic_rt_i, MS1_xic_ints_i, left_size, right_size, max_rti_1, MS1_xic_rt_t, MS1_xic_ints_t);
        List<Double> MS2_xic_rt_t = new ArrayList<>();
        List<Double> MS2_xic_ints_t = new ArrayList<>();
        trim(MS2_xic_rt_i, MS2_xic_ints_i, left_size, right_size, max_rti_2, MS2_xic_rt_t, MS2_xic_ints_t);
        // calculate scores
        double xic_dp = dot_product(MS1_xic_ints_t, MS2_xic_ints_t, MS2_xic_ints_t.size());
        double xic_norm1 = norm(MS1_xic_ints_t);
        double xic_norm2 = norm(MS2_xic_ints_t);
        double score = xic_dp / (xic_norm1 * xic_norm2);
        if (score >= XICMatchScore) {
          // put to spectra
          double spec_int = xic_dp / xic_norm1;
          SpecEntry point = new SpecEntry(mz, rt, spec_int);
          if (spectra.containsKey(MS1_xic)) {
            spectra.get(MS1_xic).add(point);
          } else {
            ArrayList<SpecEntry> points = new ArrayList<>();
            points.add(point);
            spectra.put(MS1_xic, points);
          }
          // put to ranking
          point.cor = score;
          Precursor p = new Precursor(MS1_xic.p_mz, MS1_xic.charge, MS1_xic.peak_rt, null, null);
          p.cor = score;
          if (MS1_2_cor.containsKey(p)) {
            MS1_2_cor.get(p).add(point);
          } else {
            ArrayList<SpecEntry> points = new ArrayList<>();
            points.add(point);
            MS1_2_cor.put(p, points);
          }
          if (MS2_1_cor.containsKey(point)) {
            MS2_1_cor.get(point).add(p);
          } else {
            ArrayList<Precursor> points = new ArrayList<>();
            points.add(p);
            MS2_1_cor.put(point, points);
          }
        }
      }
    }
    for (Precursor p: MS1_2_cor.keySet()) {
      Collections.sort(MS1_2_cor.get(p));
      Collections.reverse(MS1_2_cor.get(p));
    }
    for (SpecEntry p_ion: MS2_1_cor.keySet()) {
      Collections.sort(MS2_1_cor.get(p_ion));
      Collections.reverse(MS2_1_cor.get(p_ion));
    }
  }

  // use spline interpolator interpolates x and y, store in rts_i and x_i, fixed interval size
  static void interpolate(List<Double> x, List<Double> y, List<Double> x_i, List<Double> rts_i) {
    int length = x.size();
    SplineInterpolator interp = new SplineInterpolator();
    PolynomialSplineFunction f = interp.interpolate(x.stream().mapToDouble(Double::doubleValue).toArray(), y.stream().mapToDouble(Double::doubleValue).toArray());
    double start_rt = x.get(0);
    double end_rt = x.get(length - 1);
    for (double inter = start_rt; inter <= end_rt; inter += RtInter) {
      x_i.add(f.value(inter));
      rts_i.add(inter);
    }
  }

  private void trim(List<Double> x, List<Double> y, int left_s, int right_s, int center_i, List<Double> x_i, List<Double> y_i) {
    int left_bound = center_i - left_s;
    int right_bound = center_i + right_s + 1;
    x_i.addAll(x.subList(left_bound, right_bound));
    y_i.addAll(y.subList(left_bound, right_bound));
  }

  private static double dot_product(List<Double> x, List<Double> y, int len) {
    double dot_prod = 0;
    for (int d = 0; d < len; ++ d) {
      dot_prod += x.get(d) * y.get(d);
    }
    return dot_prod;
  }

  private static double norm(List<Double> x) {
    return Math.sqrt(dot_product(x, x, x.size()));
  }

  private void generate_spec() throws IOException {
    // write spectra to mgf
    int index = 0;
    for (Map.Entry<Precursor, ArrayList<SpecEntry>> spectrum : spectra.entrySet()) {
      Precursor p = spectrum.getKey();
      ArrayList<SpecEntry> xic = spectrum.getValue();
      xic.sort(Comparator.comparingDouble(specEntry -> specEntry.mz));
      info.write("BEGIN IONS\r\n");
      info.write("PEPMASS=");
      info.write(Double.toString(p.p_mz));
      info.newLine();
      info.write("CHARGE=");
      info.write(Integer.toString((int)p.charge));
      info.newLine();
      info.write("RTINSECONDS=");
      info.write(Double.toString(p.peak_rt * 60));
      info.newLine();
      info.write("TITLE=HELA-DIA-Scaffold-DIA-Trial-2-");
      info.write(Integer.toString(index));
      info.newLine();
      ArrayList<SpecEntry> points = MS1_2_cor.get(p);
      for (SpecEntry point: xic) {
        int MS1_ranking = points.indexOf(point);
        ArrayList<Precursor> ps = MS2_1_cor.get(point);
        int MS2_ranking = ps.indexOf(p);
        if (MS1_ranking < 300 && MS2_ranking < 25) {
          info.write(Double.toString(point.mz));
          info.write("\t");
          info.write(Double.toString(point.spec_int));
          info.newLine();
        }
      }
      info.write("END IONS\r\n");
      ++ index;
    }
  }
}
