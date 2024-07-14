package edu.uw.waterlooms.peptideMatch;

import java.io.Serializable;
import java.util.Arrays;

public class MSTwoTrail implements Comparable<MSTwoTrail>, Serializable {

    //mzApex	rtApex	mzs	rts	ints	peaksSum	peakArea
    private double mzApex;

    public double getRtApex() {
        return rtApex;
    }

    private double rtApex;

    public double[] getMzs() {
        return mzs;
    }


    private double intApex;

    private double[] mzs;

    public double[] getRts() {
        return rts;
    }

    private double[] rts;

    public double[] getInts() {
        return ints;
    }

    private double[] ints;
    private double peaksSum;

    public double getPeaksSum() {
        return peaksSum;
    }

    public double getPeakArea() {
        return peakArea;
    }

    private double peakArea;

    public int getPeakAreaHalfRank() {
        return peakAreaHalfRank;
    }

    public void setPeakAreaHalfRank(int peakAreaHalfRank) {
        this.peakAreaHalfRank = peakAreaHalfRank;
    }

    private int peakAreaHalfRank;

//    private List<MatchPeptdide> listMatchPeptide;

    public boolean isbSelected() {
        return (bSelected>0);
    }
    public boolean isbSelectedALL() {
        return (bSelected==1);
    }
    public boolean isbSelectedBY() {//用于检查max时，不重复检查by对应一个peak峰。
        return (bSelected==2);
    }

    public void setbSelected(int bSelected) {
        this.bSelected = bSelected;
    }

    private int bSelected = 0;

//    public void addMatchPeptide(MatchPeptdide matchPep)
//    {
//        if( listMatchPeptide == null) {
//            listMatchPeptide = new ArrayList<MatchPeptdide>();
//        }
//        listMatchPeptide.add(matchPep);
//    }
    //for data from raw data
    public MSTwoTrail(double dmz,
                      double drt,
                  double dint) {
    mzApex = dmz;
    rtApex = drt;
    peaksSum = dint;
    peakArea = dint;


    //string to stream to double[]

    mzs = new double[1];
    mzs[0] = dmz;

    rts = new double[1];
    rts[1] = drt;

    ints = new double[1];
    ints[1] = dint;

        bSelected = 0;
}

    public MSTwoTrail(double parseDouble, double parseDouble1, double[] s, double[] s1, double[] s2,
                      double parseDouble2, double parseDouble3) {
        mzApex = parseDouble;
        rtApex = parseDouble1;
        peaksSum = parseDouble2;
        peakArea = parseDouble3;


        //string to stream to double[]

        mzs = s;

        rts = s1;

        ints = s2;

        bSelected = 0;



//        mzs= Arrays.stream(s.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();
//        rts= Arrays.stream(s1.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();
//        ints= Arrays.stream(s2.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();

    }

    public MSTwoTrail(double parseDouble, double parseDouble1, String s, String s1, String s2,
                      double parseDouble2, double parseDouble3) {
        mzApex = parseDouble;
        rtApex = parseDouble1;
        peaksSum = parseDouble2;
        peakArea = parseDouble3;


        //string to stream to double[]
        String[] sp = s.split(",");
        int isp=0;
        mzs = new double[sp.length];
        for(String sp1:sp)
        {
            mzs[isp++] = Double.parseDouble(sp1);
        }
        rts = new double[sp.length];
        String[] s1p = s1.split(",");
         isp=0;
        for(String s1p1:s1p)
        {
            rts[isp++] = Double.parseDouble(s1p1);
        }
        ints = new double[sp.length];
        String[] s2p = s2.split(",");
        isp=0;
        for(String s2p1:s2p)
        {
            ints[isp++] = Double.parseDouble(s2p1);
        }
        bSelected = 0;



//        mzs= Arrays.stream(s.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();
//        rts= Arrays.stream(s1.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();
//        ints= Arrays.stream(s2.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();

    }


    public MSTwoTrail(double dmz, double drt, double dint,String s, String s1, String s2,
                      double parseDouble2, double parseDouble3) {
        mzApex = dmz;
        rtApex = drt;
        intApex = dint;
        peaksSum = parseDouble2;
        peakArea = parseDouble3;


        //string to stream to double[]
        String[] sp = s.split(",");
        int isp=0;
        mzs = new double[sp.length];
        for(String sp1:sp)
        {
            mzs[isp++] = Double.parseDouble(sp1);
        }
        rts = new double[sp.length];
        String[] s1p = s1.split(",");
        isp=0;
        for(String s1p1:s1p)
        {
            rts[isp++] = Double.parseDouble(s1p1);
        }
        ints = new double[sp.length];
        String[] s2p = s2.split(",");
        isp=0;
        for(String s2p1:s2p)
        {
            ints[isp++] = Double.parseDouble(s2p1);
        }


        bSelected = 0;

//        mzs= Arrays.stream(s.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();
//        rts= Arrays.stream(s1.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();
//        ints= Arrays.stream(s2.split(",")).flatMapToDouble(num
//                -> DoubleStream.of(Double.parseDouble(num))).toArray();

    }

    public double FindMaxIntensity() {
//        return Arrays.stream(ints).max().getAsDouble();
        int index = Utils.getMaxPos(ints);
        return ints[index];
    }


    @Override
    public int compareTo(MSTwoTrail o) {
        return Double.compare(mzApex,o.mzApex);
    }

    public double getMzApex() {
        return mzApex;
    }

    @Override
    public String toString() {
        return /*"MSTwoTrail{" +*/
                "mzApex=" + mzApex +
                ", rtApex=" + rtApex +
                ", mzs=" + Arrays.toString(mzs) +
                ", rts=" + Arrays.toString(rts) +
                ", ints=" + Arrays.toString(ints) +
                ", peaksSum=" + peaksSum +
                ", peakArea=" + peakArea +
                ", peakAreaHalfRank=" + peakAreaHalfRank ;//+
//                '}';
    }

    public MSTwoTrail(double parseDouble, double parseDouble1, double[] s, double[] s1, double[] s2,
                      double parseDouble2, double parseDouble3, int parseDouble4) {
        mzApex = parseDouble;
        rtApex = parseDouble1;
        mzs = s;
        rts = s1;
        ints = s2;
        peaksSum = parseDouble2;
        peakArea = parseDouble3;
        peakAreaHalfRank = parseDouble4;
        bSelected = 0;
    }

    public String getCurIntensity(double dms2rt) {
        for (int i =0;i< rts.length;i++)
        {
            if(rts[i]==dms2rt) {
                return String.valueOf(ints[i]);
            }
        }
        return "0";
    }

    public boolean isbSelectedMax() {
        return (bSelected==3);

    }
}
