package edu.uw.waterlooms.peptideMatch;

import org.json.JSONWriter;

import java.util.ArrayList;
import java.util.Arrays;

public class MSOneTrail implements Comparable<MSOneTrail>{
    public static double threshold = 0.005; // mz error allowance



    public long getId() {
        return id;
    }

    //id	mz	rt	z	isotope_num	intensity_shape_score	isotope_distribution_score	intensity_window_evg	intensity_area_percentage	rt_start	rt_end	scan_num	intensity_sum	svr_score	quality_score
    private long id;

    public double getMz() {
        return mz;
    }

    public double getMass()
    {
        return Utils.MzToMass(mz,z);
    }
    private double mz;


    public double getStratRt() {
        return arrMS1Rts[0];
    }
    public double getEndRt() {
        return arrMS1Rts[arrMS1Rts.length-1];
    }

    public boolean isIzCombine() {
        return izCombine;
    }

    public void setIzCombine(boolean izCombine) {
        this.izCombine = izCombine;
    }

    private boolean izCombine;


    public double getRt() {
        return rt;
    }

    private double rt;

    public int getZ() {
        return z;
    }

    private int z;
    private int isotope_num;
    private double intensity_shape_score;
    private double isotope_distribution_score;
    private double  intensity_window_evg;
    private double intensity_area_percentage;
    private double rt_start;
    private double rt_end;

    public long getScan_num() {
        return scan_num;
    }

    private long scan_num;

    private double intensity_sum;
    private double svr_score;

    private double quantification_peaks_sum;

    public double getQuantification_peaks_area() {
        return quantification_peaks_area;
    }

    private double quantification_peaks_area;
    private String mzs;
    private String rts;
    private String ints;

    String strGeneratePeptide;

    public boolean isbAddPrecursor() {
        return bAddPrecursor;
    }

    boolean bAddPrecursor;
    public double[] getArrMS1Mzs() {
        return arrMS1Mzs;
    }

    public void setArrMS1Mzs(ArrayList<Double> listMS1Mzs) {
        this.arrMS1Mzs = listMS1Mzs.stream().mapToDouble(Double::doubleValue).toArray();
    }

    private double[] arrMS1Mzs;

    public double[] getArrMS1Rts() {
        return arrMS1Rts;
    }

    public void setArrMS1Rts(ArrayList<Double> listMS1Rts) {
        this.arrMS1Rts = listMS1Rts.stream().mapToDouble(Double::doubleValue).toArray();
    }

    private double[] arrMS1Rts;

    public double[] getArrMS1Ins() {
        return arrMS1Ins;
    }

    public void setArrMS1Ins(ArrayList<Double> listMS1Ins) {
        this.arrMS1Ins = listMS1Ins.stream().mapToDouble(Double::doubleValue).toArray();
    }

    private double[] arrMS1Ins;

    public double getQuality_score() {
        return quality_score;
    }

    private double quality_score;

    private double dmaxmz;

    private double dminmz;

    public void setdPeakAreaLocalRankInTimeSpan(double dPeakAreaLocalRankInTimeSpan) {
        this.dPeakAreaLocalRankInTimeSpan = dPeakAreaLocalRankInTimeSpan;
    }

    public double getdPeakAreaLocalRankInTimeSpan() {
        return dPeakAreaLocalRankInTimeSpan;
    }

    private double dPeakAreaLocalRankInTimeSpan;


    public void quantifyPeaks() {

        double peaks_sum = arrMS1Ins[0];
        double peaks_area = 0;
        for (int j = 1; j < arrMS1Ins.length; j++) {
            peaks_sum += arrMS1Ins[j];
            peaks_area += (arrMS1Ins[j] + arrMS1Ins[j-1]) * (arrMS1Rts[j] - arrMS1Rts[j-1]) / 2;
        }
        //若只有一个峰，peakarea则负值为该峰的丰度*时间差/2，即后一点为0
        if (arrMS1Ins.length==1)
        {
            peaks_area = arrMS1Ins[0]*Utils.MS1RTSpan/2;
        }
        quantification_peaks_sum = peaks_sum;
        quantification_peaks_area = peaks_area;
    }

    public MSOneTrail(long parseDouble, double parseDouble1, double parseDouble2, int parseDouble3, int parseDouble4,
                      double parseDouble5, double parseDouble6, double parseDouble7, double parseDouble8, double parseDouble9,
                      double parseDouble10, long parseDouble11, double parseDouble12, double parseDouble13, double parseDouble14)
    {
        id = parseDouble;
        mz = parseDouble1;
        rt = parseDouble2;
        z = parseDouble3;
        isotope_num = parseDouble4;
        intensity_shape_score = parseDouble5;
        isotope_distribution_score = parseDouble6;
        intensity_window_evg = parseDouble7;
        intensity_area_percentage = parseDouble8;
        rt_start = parseDouble9;
        rt_end = parseDouble10;
        scan_num = parseDouble11;
        intensity_sum = parseDouble12;
        svr_score = parseDouble13;
        quality_score = parseDouble14;
        bAddPrecursor = false;
    }

    public MSOneTrail(long lid, double dmz, double drt,
                      int iz, int iisotope_num, long lscan_num,
                      double dquantification_peaks_area) {
        id = lid;
        mz = dmz;
        rt = drt;
        z = iz;
        isotope_num = iisotope_num;
        intensity_shape_score = 0;
        isotope_distribution_score = 0;
        intensity_area_percentage = 0;
        scan_num = lscan_num;
        quantification_peaks_sum =  dquantification_peaks_area;
        quantification_peaks_area = dquantification_peaks_area;
        svr_score = 0;
        quality_score = 0;

        arrMS1Mzs = new double[1];
        arrMS1Mzs[0] = dmz;

        arrMS1Rts = new double[1];
        arrMS1Rts[0] = drt;

        arrMS1Ins = new double[1];
        arrMS1Ins[0] = dquantification_peaks_area;
        bAddPrecursor = false;

    }

    public MSOneTrail(long lid, double dmz, double drt,
                      int iz, int iisotope_num, double dintensity_shape_score,
                      double disotope_distribution_score, double dintensity_area_percentage, long lscan_num,
                      double dquantification_peaks_sum, double dquantification_peaks_area, double dsvr_score,
                      double dquality_score, String smz, String srt, String sints) {
        id = lid;
        mz = dmz;
        rt = drt;
        z = iz;
        isotope_num = iisotope_num;
        intensity_shape_score = dintensity_shape_score;
        isotope_distribution_score = disotope_distribution_score;
        intensity_area_percentage = dintensity_area_percentage;
        scan_num = lscan_num;
        quantification_peaks_sum =  dquantification_peaks_sum;
        quantification_peaks_area = dquantification_peaks_area;
        svr_score = dsvr_score;
        quality_score = dquality_score;
        if(smz.charAt(0)=='[')
        {
            mzs = smz.substring(1,smz.length()-1);
            rts	= srt.substring(1,srt.length()-1);
            ints = sints.substring(1,sints.length()-1);
        }else
        {
            mzs = smz;
            rts	= srt;
            ints = sints;
        }



        String[] sp = mzs.split(",");
        int isp=0;
        arrMS1Mzs = new double[sp.length];
        for(String sp1:sp)
        {
            arrMS1Mzs[isp++] = Double.parseDouble(sp1);
            if(isp==1)
            {
                dmaxmz = arrMS1Mzs[isp-1];
                dminmz = arrMS1Mzs[isp-1];

            }else {
                if(arrMS1Mzs[isp-1]>dmaxmz)
                    dmaxmz = arrMS1Mzs[isp-1];
                if(arrMS1Mzs[isp-1]<dminmz)
                    dminmz = arrMS1Mzs[isp-1];
            }

        }
        arrMS1Rts = new double[sp.length];
        String[] s1p = rts.split(",");
         isp=0;
        for(String s1p1:s1p)
        {
            arrMS1Rts[isp++] = Double.parseDouble(s1p1);
        }
        arrMS1Ins = new double[sp.length];
        String[] s2p = ints.split(",");
        isp=0;
        for(String s2p1:s2p)
        {
            arrMS1Ins[isp++] = Double.parseDouble(s2p1);
        }
        bAddPrecursor = false;
    }

    //判断是不是在该fragment的trail中。条件一是在rt范围内，条件二是mz大于ppm/2的最小值且小于ppm/2的最大值
    public boolean isRTMZinRTSMZS(double tmprt,double tmpmz,int iz)
    {
        boolean isIn = true;
        if(iz!=z)
        {
            isIn = false;
            return isIn;
        }
        if(tmprt<arrMS1Rts[0] || tmprt>arrMS1Rts[arrMS1Rts.length-1])
        {
            isIn = false;
            return isIn;
        }

//        if (tmpmz<mz*(1-Utils.thresholdPPM/2) || tmpmz>mz*(1+Utils.thresholdPPM/2))
        if (tmpmz<dminmz || tmpmz>dmaxmz)
        {
            isIn = false;
            return isIn;
        }
        return isIn;
    }
    @Override
    public int compareTo(MSOneTrail o) {
        return Double.compare(mz, o.getMz());

    }


    @Override
    public String toString() {
        return   ""+ id +
                "\t" + mz +
                "\t" + rt +
                "\t" + z +
                "\t" + isotope_num +
                "\t" + intensity_shape_score +
                "\t" + isotope_distribution_score +
                "\t" + intensity_area_percentage +
                "\t" + scan_num +
                "\t" + quantification_peaks_sum +
                "\t" + quantification_peaks_area +
                "\t" + svr_score +
                "\t" + quality_score +
                "\t" + JSONWriter.valueToString(arrMS1Mzs) +
                "\t" + JSONWriter.valueToString(arrMS1Rts) +
                "\t" + JSONWriter.valueToString(arrMS1Ins) +
                "\t" + bAddPrecursor +
                "\n";
    }

    public void setQualityScore(int i) {
        quality_score = (i + isotope_num)/2.0;
    }
}
