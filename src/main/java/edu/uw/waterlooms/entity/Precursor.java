package edu.uw.waterlooms.entity;

import java.util.List;

public class Precursor implements Comparable<Precursor> {
  public double p_mz;
  public double charge;
  public double peak_rt;
  public List<Double> xic_rt;
  public List<Double> xic_ints;
  public double cor = 0;

  public Precursor(double d1, double d2, double d3, List<Double> d4, List<Double> d5) {
    p_mz = d1;
    charge = d2;
    peak_rt = d3;
    xic_rt = d4;
    xic_ints = d5;
  }

  public boolean equals(Object obj) {
    if (obj == null) return false;
    if (obj == this) return true;
    if (!(obj instanceof Precursor)) return false;
    Precursor o = (Precursor) obj;
    return o.p_mz == this.p_mz && o.charge == this.charge && o.peak_rt == this.peak_rt;
  }

  @Override
  public int compareTo(Precursor other) {
    return Double.compare(cor, other.cor);
  }

  @Override
  public int hashCode() {
    return (Double.valueOf(this.p_mz).hashCode()
        + Double.valueOf(this.charge).hashCode()
        + Double.valueOf(this.peak_rt).hashCode());
  }
}
