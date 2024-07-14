package edu.uw.waterlooms.peptideMatch;


//import jdk.nashorn.internal.ir.debug.JSONWriter;
import org.apache.commons.math3.util.FastMath;
import org.json.JSONWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

import static java.lang.Math.abs;
import static java.lang.Math.log10;

public class DbMatch implements Comparable<DbMatch>{
    public MatchedInterval matchedRts; // stores the interval where peptide most likely occurs
//    public Double score;
    public String id; // id of protein
    public String composition; // peptide composition
//    static List<Integer> bpos = new ArrayList<>();
//    static List<Integer> ypos = new ArrayList<>();
//
//    static List<Integer> bipos = new ArrayList<>();
//    static List<Integer> yipos = new ArrayList<>();



    static int debugOutput =1;
//    static List<Double> bMatchDistance = new ArrayList<>();
//    static List<Double> yMatchDistance = new ArrayList<>();





    public DbMatch(String id, String composition, MatchedInterval matchedRts) {
        this.composition = composition;
        this.matchedRts = matchedRts;
        this.id = id;
    }

    @Override
    public int compareTo(DbMatch other) {
        Double score = this.matchedRts.score;
        return score.compareTo(other.matchedRts.score);
    }

    public static double EmpiricalScore(double intensity) {
        // so score is >=0
        if (intensity < 1) {
            return 0;
        }
        return log10(intensity);

//        return 1;
    }

    // returns if peak is already in rtSet
    public static boolean checkDup(ArrayList<Peak> rtSet, Peak peak, Integer index) {
        if (rtSet.size() == 1) { // index would be -1
            if (rtSet.get(-(index+1)).checkSame(peak)) {
                return true;
            }
        }

        if (index > 0 && index < rtSet.size()) { // check rtSet[i], rtSet[i-1]
            if (rtSet.get(index-1).checkSame(peak)) { // check the one before
                return true;
            }
            if (rtSet.get(-(index+1)).checkSame(peak)) {
                return true;
            }
        }
        return false;
    }

//SELECTED
    public  static  MatchPeptideWithMultiCharge PepMatchIonCountCrossRTMultiChargeFragmentLoss_top6PredictIntensity(Peptide pep, List<MSTwoTrail> lMSTwotrails, double maxPeakArea,
                                                                                                                    Enums.ScoreMode scoreMode, BufferedWriter writer, int icount,
                                                                                                                    MSOneTrail msOneTrail/*double ms1rt, double ms1Mass, long ms1id,double ms1qualityScore*/,
                                                                                                                    MSTwoTrail[] bIonMatchTrail, MSTwoTrail[] yIonMatchTrail, int[] arrbloss, int[] arryloss,
                                                                                                                    double[][] arrIntensity, MS2Intensity ms2Int) {
    List<Integer> bpos = null;
    List<Integer> ypos = null;
    List<Integer> bNormalpos  = Collections.synchronizedList(new ArrayList<>());
    List<Integer> yNormalpos  = Collections.synchronizedList(new ArrayList<>());
    List<Integer> bH2OLosspos = Collections.synchronizedList(new ArrayList<>());
    List<Integer> yH2OLosspos = Collections.synchronizedList(new ArrayList<>());
    List<Integer> bNH3Losspos = Collections.synchronizedList(new ArrayList<>());
    List<Integer> yNH3Losspos = Collections.synchronizedList(new ArrayList<>());

    List<Integer> bipos = Collections.synchronizedList(new ArrayList<>());
    List<Integer> yipos = Collections.synchronizedList(new ArrayList<>());

//        List<Double> bMatchDistance =  Collections.synchronizedList(new ArrayList<>());
//        List<Double> yMatchDistance = Collections.synchronizedList(new ArrayList<>());

    double dMatch = 1.0;
    double[] dlogBionPeakAreaSum = new double[Utils.iMS2Charge];
    double[] dlogYionPeakAreaSum = new double[Utils.iMS2Charge];
    double[] dBionMassErrorSum = new double[Utils.iMS2Charge];
    double[] dYionMassErrorSum = new double[Utils.iMS2Charge];
    double[] dBionRetentionTimeErrorSum = new double[Utils.iMS2Charge];
    double[] dYionRetentionTimeErrorSum = new double[Utils.iMS2Charge];
    double[] dBionPeakHalfSum = new double[Utils.iMS2Charge];
    double[] dYionPeakHalfSum = new double[Utils.iMS2Charge];
    int[] iBConsective = new int[Utils.iMS2Charge];
    int[] iYConsective = new int[Utils.iMS2Charge];
    String strb = "";
    String stry = "";
    int[] bSum = new int[Utils.iMS2Charge];
    int[] ySum = new int[Utils.iMS2Charge];
    int bySum = 0;
    int iBYComplementaryTrailCount = 0;

    double sumXX = 0d;
    double sumYY = 0d;
    double sumXY = 0d;

    double sumXXMatchPredict = 0d;
    double sumYYMatchPredict = 0d;
    double sumXYMatchPredict = 0d;

    double db60sumXX = 0d;
    double db60sumYY = 0d;
    double db60sumXY = 0d;

    double dbTop6sumXX = 0d;
    double dbTop6sumYY = 0d;
    double dbTop6sumXY = 0d;

    double dyTop6sumXX = 0d;
    double dyTop6sumYY = 0d;
    double dyTop6sumXY = 0d;

    double dcharget1Top6sumXX = 0d;
    double dcharget1Top6sumYY = 0d;
    double dcharget1Top6sumXY = 0d;

    double dallTop6sumXX = 0d;
    double dallTop6sumYY = 0d;
    double dallTop6sumXY = 0d;

    int itopall= 0;
    int ibc1Top6= 0;
    int iyc1Top6= 0;
    int iallc1Top6= 0;

    TreeMap<Integer, LinkedList<MSTwoTrail>>  mapBionTrail = null;
    TreeMap<Integer, LinkedList<MSTwoTrail>>  mapYionTrail = null;

//        TreeMap<Integer, MSTwoTrail>  mapCurRTBionMS2Trail = null;
//        TreeMap<Integer, MSTwoTrail>  mapCurRTYionMS2Trail = null;

    double[] dBionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];
    double[] dYionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];

    int ibAllMatchedWithloss = 0;
    int ibMatchNormalAndH2Oloss = 0;
    int ibMatchNormalAndNH3loss = 0;
    int iyAllMatchedWithloss = 0;
    int iyMatchNormalAndH2Oloss = 0;
    int iyMatchNormalAndNH3loss = 0;

    if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances || Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
    {
        dMatch = 0.0;//add all abundance
    }
    int len = pep.b_ions.length;
    MatchPeptideWithMultiCharge matchPep = null;

//        int iCurrentPos = 0;
    for(int iCharge=1;iCharge<=Utils.iMS2Charge;iCharge++) {
//            bpos.clear();
//            ypos.clear();
        bNormalpos.clear();
        yNormalpos.clear();
        bH2OLosspos.clear();
        yH2OLosspos.clear();
        bNH3Losspos.clear();
        yNH3Losspos.clear();

        bipos.clear();
        yipos.clear();


        int ibposFront = 0;
        int ibposBack = 0;
        int iyposFront = 0;
        int iyposBack = 0;

        int ibnearpos = 0;
        int iynearpos = 0;

        int bMatchSize;
        int yMatchSize;

        double dbsumXX = 0d;
        double dbsumYY = 0d;
        double dbsumXY = 0d;

        double dysumXX = 0d;
        double dysumYY = 0d;
        double dysumXY = 0d;




        for (int i = 0; i < len -1; ++i) {// b y  ion less than len of peptides

            //发现fragment loos有两种策略，一是将丢失的离子看成是一个正常的可匹配的离子峰，二是将其作为已匹配的正常离子峰的附加属性（若匹配的正常的离子再检查他有没有附加的丢失水和胺的离子峰并统计其个数）
            //下面的代码是策略一
            int ibCountFragmentTypeMatched = 0;
            int iyCountFragmentTypeMatched = 0;

            for (int iFragmentType = 0; iFragmentType < Utils.iMS2FragmentTypeCount; iFragmentType++) {

                double b_ionsMZ = (pep.b_ions[i] + AminoAcid.ProtonMass * (iCharge - 1)) / iCharge;//由于b离子原来加过一个离子，对后面的离子只需要加入一个质子质量 即可 mass to mz 该公式化简是 mass/z + AminoAcid.ProtonMass;
                double y_ionsMZ = (pep.y_ions[len - i - 1] + AminoAcid.ProtonMass * (iCharge - 1)) / iCharge;//由于y离子原来加过一个离子，对后面的离子只需要加入一个质子质量 即可
                if(Utils.iMS2FragmentTypeCount==3) {//考虑loss
                    //从丢失的 大小开始，越大越在前面，mz越小，这样search时候先搜小的然后搜大的
                    if (iFragmentType == 0) {
                        bpos = bH2OLosspos;
                        ypos = yH2OLosspos;
                        b_ionsMZ = b_ionsMZ - (Utils.H2OMass) / iCharge;
                        y_ionsMZ = y_ionsMZ - (Utils.H2OMass) / iCharge;
                    } else if (iFragmentType == 1) {
                        bpos = bNH3Losspos;
                        ypos = yNH3Losspos;
                        b_ionsMZ = b_ionsMZ - (Utils.NH3Mass) / iCharge;
                        y_ionsMZ = y_ionsMZ - (Utils.NH3Mass) / iCharge;
                    } else {
                        bpos = bNormalpos;
                        ypos = yNormalpos;
                    }
                }else if(Utils.iMS2FragmentTypeCount==1)//不考虑loss
                {
                    bpos = bNormalpos;
                    ypos = yNormalpos;
                }
                bMatchSize = bpos.size();
                yMatchSize = ypos.size();

//                double minBIonDistance = Utils.thresholdPPM;
//                double minYIonDistance = Utils.thresholdPPM;

                double[][] arr_lftIsotope = new double[Utils.iMS2Charge][2];

                double[][] arr_rightIsotope = new double[Utils.iMS2Charge][2];
                for (int ic = 0; ic < Utils.iMS2Charge; ic++) {
                    arr_lftIsotope[ic][0] = b_ionsMZ * (1 - Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                    arr_lftIsotope[ic][1] = b_ionsMZ * (1 + Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                    arr_rightIsotope[ic][0] = b_ionsMZ * (1 - Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                    arr_rightIsotope[ic][1] = b_ionsMZ * (1 + Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                }

                double bIonMinusMass = b_ionsMZ * (1 - Utils.thresholdPPM);
                double bIonMaxMass = b_ionsMZ * (1 + Utils.thresholdPPM);

                int iBionCurCount = 0;
                int iYionCurCount = 0;


                //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
                if ((bIonMaxMass > (lMSTwotrails.get(0).getMzApex()))
                        && (bIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size() - 1).getMzApex()))) {
                    ibposFront = binarySearch0(lMSTwotrails, ibposFront, lMSTwotrails.size(), bIonMinusMass);//全新查找
                    ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposFront);

                    ///更新bpos里面的内容
//                    iBionCurCount = getiIonCurCount(lMSTwotrails, ibnearpos, bIonMaxMass, bpos);
//                    iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass,b_ionsPreIsotopeMaxMZ,b_ionsPreIsotopeMinMZ
//                            ,arr_rightIsotope,iCharge, bpos);
                    iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass, arr_lftIsotope
                            , arr_rightIsotope, iCharge, bpos);

                    if (bpos.size() > 0) {
                        ibposFront = bpos.get(bpos.size() - 1);
                    }
                    if (ibposFront < 0) {
                        ibposFront = ibnearpos;
                    }
                }


//                double y_ionsPreIsotopeMinMZ = y_ionsMZ * (1 -  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / iCharge;
//                double y_ionsPreIsotopeMaxMZ = y_ionsMZ * (1 +  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / iCharge;


                for (int ic = 0; ic < Utils.iMS2Charge; ic++) {
                    arr_lftIsotope[ic][0] = y_ionsMZ * (1 - Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                    arr_lftIsotope[ic][1] = y_ionsMZ * (1 + Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                    arr_rightIsotope[ic][0] = y_ionsMZ * (1 - Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                    arr_rightIsotope[ic][1] = y_ionsMZ * (1 + Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                }


                double yIonMinusMass = y_ionsMZ * (1 - Utils.thresholdPPM);
                double yIonMaxMass = y_ionsMZ * (1 + Utils.thresholdPPM);


                //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
                if ((yIonMaxMass > (lMSTwotrails.get(0).getMzApex()))
                        && (yIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size() - 1).getMzApex()))) {
                    //根据前面比较的结果，search 对应的y离子

                    iyposFront = binarySearch0(lMSTwotrails, iyposFront, lMSTwotrails.size(), yIonMinusMass);
                    iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);

//                    iYionCurCount = getiIonCurCount(lMSTwotrails, iynearpos, yIonMaxMass, ypos);
//                    iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass,y_ionsPreIsotopeMaxMZ,y_ionsPreIsotopeMinMZ
//                            ,arr_rightIsotope,iCharge, ypos);
                    iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass, arr_lftIsotope
                            , arr_rightIsotope, iCharge, ypos);


                    if (ypos.size() > 0) {
                        iyposFront = ypos.get(ypos.size() - 1);
                    }
                    if (iyposFront < 0) {
                        iyposFront = iynearpos;
                    }
                }


                //判断单个bion with ms2 trail匹配是否成功
                //判断连续ion
                //大于某阈值后的连续离子，大于ms2的中平均数？中位数？大于peptide的质量的一半？
                //根据ms2 trail数据的散列度设置阈值？

                //用于记录哪个离子峰匹配，包含：loss类型、电荷、离子峰编号
                //framentType+icharge+iNumber,没有loss排前面
                int currentPosWithMatch = (Utils.iMS2FragmentTypeCount-iFragmentType-1) * (Utils.iMS2Charge * len) +(iCharge - 1) * len + i  ;
                if (bpos.size() > bMatchSize) {
                    ibCountFragmentTypeMatched += (iFragmentType + 1);//建立特征--不考虑跨rt
//                        if (mapCurRTBionMS2Trail == null) {
//                            mapCurRTBionMS2Trail = new TreeMap<>();
//                        }
//                        mapCurRTBionMS2Trail.put(currentPosWithMatch + 1, lMSTwotrails.get(ibposFront));//其他地方没有用到
                    if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        bIonMatchTrail[currentPosWithMatch] = lMSTwotrails.get(ibposFront);
                    } else//匹配多个MS2RT
                    {
                        bIonMatchTrail[currentPosWithMatch] = CompateTwoTrail(lMSTwotrails.get(ibposFront), bIonMatchTrail[currentPosWithMatch], b_ionsMZ);
                    }
                } else {
                    if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        bIonMatchTrail[currentPosWithMatch] = null;
                    }
                }

                if (bIonMatchTrail[currentPosWithMatch] != null)//bpos.size()>bMatchSize)
                {
                    if ((bipos.size() > 0) && (bipos.get(bipos.size() - 1) == currentPosWithMatch - 1)) {
                        iBConsective[iCharge - 1]++;
                    }

                    bipos.add(currentPosWithMatch);

                    if (mapBionTrail == null) {
                        mapBionTrail = new TreeMap<>();
                    }
                    LinkedList<MSTwoTrail> listMS2Trail = mapBionTrail.get(currentPosWithMatch);
                    if (listMS2Trail == null) {
                        listMS2Trail = new LinkedList<MSTwoTrail>();
                    }
                    listMS2Trail.add(bIonMatchTrail[currentPosWithMatch]);
                    mapBionTrail.put(currentPosWithMatch, listMS2Trail);

                    dMatch = getdMatchIonScore(bIonMatchTrail[currentPosWithMatch], bIonMaxMass, dMatch, maxPeakArea, iBionCurCount);

                    dMatch = getdMatchMultiIonScore(lMSTwotrails, bipos, bpos, dMatch);

                    //用于输出中间值

                    dlogBionPeakAreaSum[iCharge - 1] += Math.log10(bIonMatchTrail[currentPosWithMatch].getPeakArea());
                    dBionMassErrorSum[iCharge - 1] += 1.0 - Math.pow((((bIonMatchTrail[currentPosWithMatch].getMzApex() - b_ionsMZ) / (bIonMatchTrail[currentPosWithMatch].getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                    dBionPeakHalfSum[iCharge - 1] += Math.max(Math.log10(50.0 / bIonMatchTrail[currentPosWithMatch].getPeakAreaHalfRank()), 0.0);
                    dBionRetentionTimeErrorSum[iCharge - 1] += 1 - Math.pow((bIonMatchTrail[currentPosWithMatch].getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
                    if (arrIntensity != null) {
//                        System.out.println("B arrIntensity "+pep.composition+"  "+msOneTrail.getZ());

                        double dxy = bIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge - 1];
                        double dyy = arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];
                        double dxx = bIonMatchTrail[currentPosWithMatch].getPeakArea() * bIonMatchTrail[currentPosWithMatch].getPeakArea();
                        sumXY += dxy;
                        sumYY += dyy;
                        sumXX += dxx;

                        dbsumXY += dxy;
                        dbsumYY += dyy;
                        dbsumXX += dxx;

                        if (iCharge == 1 && i < len * 0.6) {
                            db60sumXY += dxy;
                            db60sumYY += dyy;
                            db60sumXX += dxx;

                        }
                        //若在预测丰度的top6里的b离子
//                        if(iCharge == 1 && i== ms2Int.arrbTop6IntensityPosCharge1[ibtop]  )
                        if(ms2Int!=null && ms2Int.arrallTop6IntensityPos!=null)
                        {
                            if(itopall<Utils.iNumPredictIntensityTopN && ms2Int.arrallTop6IntensityCount[iCharge-1]>0)
                            {
                                if(ms2Int.arrallTop6IntensityPos[itopall]==i)
                                {
                                    dallTop6sumXY += dxy;
                                    dallTop6sumYY += dyy;
                                    dallTop6sumXX += dxx;
                                    itopall++;
                                }
                            }

                            if(iCharge ==1 && iallc1Top6<Utils.iNumPredictIntensityTopN &&  ms2Int.arrallTop6IntensityCharge1>iallc1Top6)
                            {
                                if(ms2Int.arrallTop6IntensityCharge1Pos[iallc1Top6]==i)
                                {
                                    dcharget1Top6sumXY += dxy;
                                    dcharget1Top6sumYY += dyy;
                                    dcharget1Top6sumXX += dxx;
                                    iallc1Top6++;
                                }
                            }


                            if(iCharge ==1 && ibc1Top6<Utils.iNumPredictIntensityTopN &&  ms2Int.arrbTop6IntensityPosCharge1[ibc1Top6]==i)
                            {
                                dbTop6sumXY += dxy;
                                dbTop6sumYY += dyy;
                                dbTop6sumXX += dxx;
                                ibc1Top6++;
                            }
                        }


                        if (arrIntensity[i][iCharge - 1] > 10 * Double.MIN_VALUE) {
                            sumXYMatchPredict += dxy;
                            sumYYMatchPredict += dyy;
                            sumXXMatchPredict += dxx;

                        }

                        //debug
                        if (msOneTrail.getId()==157058 && pep.composition.equals("VVAGVAAALAHK") && ms2Int!=null
                                && ms2Int.arrallTop6IntensityPos!=null){
                            if(Utils.iDebugoutput > 0)
                            {
                                try {
                                    Utils.bufferedWriterDebugInfo.write("rt2:"+ bIonMatchTrail!=null && bIonMatchTrail.length>0 && bIonMatchTrail[currentPosWithMatch]!=null? ""+ bIonMatchTrail[currentPosWithMatch].getRtApex() :"0"+"\n");
                                    Utils.bufferedWriterDebugInfo.write("bi:"+ i+"\n");
                                    Utils.bufferedWriterDebugInfo.write("bIonMatchTrail[currentPosWithMatch].getPeakArea() :"+bIonMatchTrail[currentPosWithMatch].getPeakArea() +"\n");
                                    Utils.bufferedWriterDebugInfo.write("arrIntensity[i][iCharge - 1]:"+arrIntensity[i][iCharge - 1]+"\n");


                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }


                            }

                        }
                    }
//                    else
//                    {
//                        if(debugOutput++<10) {
//                            System.out.println(" no arrIntensity "+pep.composition+"  "+msOneTrail.getZ());
//                        }
//                    }

                } else {

                    if (arrIntensity != null) {
                        double dyy = arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];
//                        System.out.println("   bnotin----"+i+"  "+pep.composition+"   "+ msOneTrail.getZ());

//                        System.out.println("arrIntensity:" + JSONWriter.valueToString(arrIntensity));
//                        sumXY += Utils.NoneMatchedIntensity * arrIntensity[i][iCharge -1 ];
//                        sumXX += Utils.NoneMatchedIntensity * Utils.NoneMatchedIntensity;
                        sumYY += dyy;//做cosine时，没有匹配Y值（不是y离子）也要加

                        dbsumYY += dyy;//做cosine时，没有匹配Y值（不是y离子）也要加

                        if (iCharge == 1 && i < len * 0.6) {
                            db60sumYY += dyy;//做cosine时，没有匹配Y值（不是y离子）也要加

                        }


                        if(ms2Int!=null && ms2Int.arrallTop6IntensityPos!=null)
                        {
                            if(itopall<Utils.iNumPredictIntensityTopN && ms2Int.arrallTop6IntensityCount[iCharge-1]>0)
                            {
                                if(ms2Int.arrallTop6IntensityPos[itopall]==i)
                                {
                                    dallTop6sumYY += dyy;
                                    itopall++;
                                }
                            }

                            if(iCharge ==1 && iallc1Top6<Utils.iNumPredictIntensityTopN &&  ms2Int.arrallTop6IntensityCharge1>iallc1Top6)
                            {
                                if(ms2Int.arrallTop6IntensityCharge1Pos[iallc1Top6]==i)
                                {
                                    dcharget1Top6sumYY += dyy;
                                    iallc1Top6++;
                                }
                            }


                            if(iCharge ==1 && ibc1Top6<Utils.iNumPredictIntensityTopN &&  ms2Int.arrbTop6IntensityPosCharge1[ibc1Top6]==i)
                            {
                                dbTop6sumYY += dyy;
                                ibc1Top6++;
                            }
                        }
                        //debug
                        if (msOneTrail.getId()==157058 && pep.composition.equals("VVAGVAAALAHK") && ms2Int!=null
                                && ms2Int.arrallTop6IntensityPos!=null){
                            if(Utils.iDebugoutput > 0)
                            {
                                try {
                                    Utils.bufferedWriterDebugInfo.write("bi:"+ i+"\n");
                                    Utils.bufferedWriterDebugInfo.write("arrIntensity[i][iCharge - 1]:"+arrIntensity[i][iCharge - 1]+"\n");


                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }


                            }

                        }
                    }
                }

                //用于记录哪个离子峰匹配，包含：loss类型、电荷、离子峰编号
                //framentType+icharge+iNumber,没有loss排前面
                currentPosWithMatch = (Utils.iMS2FragmentTypeCount-iFragmentType-1) * (Utils.iMS2Charge * len) + (iCharge - 1) * len + len - i - 1  ;

                if (ypos.size() > yMatchSize) {
                    iyCountFragmentTypeMatched += (iFragmentType + 1);//建立特征--不考虑跨rt

//                        if (mapCurRTYionMS2Trail == null) {
//                            mapCurRTYionMS2Trail = new TreeMap<>();
//                        }
//                        mapCurRTYionMS2Trail.put((iCharge - 1) * len + i + 1, lMSTwotrails.get(iyposFront));//其他地方没有用到
                    if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        yIonMatchTrail[currentPosWithMatch] = lMSTwotrails.get(iyposFront);
                    } else//匹配多个MS2RT
                    {
                        yIonMatchTrail[currentPosWithMatch] = CompateTwoTrail(lMSTwotrails.get(iyposFront), yIonMatchTrail[currentPosWithMatch], y_ionsMZ);
                    }
                } else {
                    if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        yIonMatchTrail[currentPosWithMatch] = null;
                    }
                }

                if (yIonMatchTrail[currentPosWithMatch] != null) {
                    if ((yipos.size() > 0) && (yipos.get(yipos.size() - 1) == (iCharge - 1) * len + len - i)) {
                        iYConsective[iCharge - 1]++;
                    }


                    yipos.add(currentPosWithMatch);

                    if (mapYionTrail == null) {
                        mapYionTrail = new TreeMap<>();
                    }
                    LinkedList<MSTwoTrail> listMS2Trail = mapYionTrail.get(currentPosWithMatch);
                    if (listMS2Trail == null) {
                        listMS2Trail = new LinkedList<MSTwoTrail>();
                    }
                    listMS2Trail.add(yIonMatchTrail[currentPosWithMatch]);
                    mapYionTrail.put(currentPosWithMatch, listMS2Trail);


                    dMatch = getdMatchIonScore(yIonMatchTrail[currentPosWithMatch], yIonMaxMass, dMatch, maxPeakArea, iYionCurCount);
                    dMatch = getdMatchMultiIonScore(lMSTwotrails, yipos, ypos, dMatch);

                    dlogYionPeakAreaSum[iCharge - 1] += Math.log10(yIonMatchTrail[currentPosWithMatch].getPeakArea());
                    dYionMassErrorSum[iCharge - 1] += 1.0 - Math.pow((((yIonMatchTrail[currentPosWithMatch].getMzApex() - y_ionsMZ) / (yIonMatchTrail[currentPosWithMatch].getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                    dYionPeakHalfSum[iCharge - 1] += Math.max(Math.log10(50.0 / yIonMatchTrail[currentPosWithMatch].getPeakAreaHalfRank()), 0.0);
                    dYionRetentionTimeErrorSum[iCharge - 1] += 1 - Math.pow((yIonMatchTrail[currentPosWithMatch].getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
                    if (arrIntensity != null) {
//                        System.out.println("Y arrIntensity "+pep.composition+"  "+msOneTrail.getZ());
                        double dxy = yIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge + 1];
                        double dyy = arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                        double dxx = yIonMatchTrail[currentPosWithMatch].getPeakArea() * yIonMatchTrail[currentPosWithMatch].getPeakArea();

                        sumXY += dxy;
                        sumYY += dyy;
                        sumXX += dxx;

                        dysumXY += dxy;
                        dysumYY += dyy;
                        dysumXX += dxx;

                        if(ms2Int!=null && ms2Int.arrallTop6IntensityPos!=null) {
                            if (itopall < Utils.iNumPredictIntensityTopN && ms2Int.arrallTop6IntensityCount[iCharge + 1] > 0) {
                                if (ms2Int.arrallTop6IntensityPos[itopall] == i) {
                                    dallTop6sumXY += dxy;
                                    dallTop6sumYY += dyy;
                                    dallTop6sumXX += dxx;

                                    itopall++;

                                }
                            }


                            if (iCharge == 1 && iallc1Top6 < Utils.iNumPredictIntensityTopN && ms2Int.arrallTop6IntensityCharge1 <= iallc1Top6) {
                                if (ms2Int.arrallTop6IntensityCharge1Pos[iallc1Top6] == i ) {
                                    dcharget1Top6sumXY += dxy;
                                    dcharget1Top6sumYY += dyy;
                                    dcharget1Top6sumXX += dxx;

                                    iallc1Top6++;

                                }
                            }


                            if (iCharge == 1 && iyc1Top6 < Utils.iNumPredictIntensityTopN && ms2Int.arryTop6IntensityPosCharge1[iyc1Top6] ==  i ) {
                                dyTop6sumXY += dxy;
                                dyTop6sumYY += dyy;
                                dyTop6sumXX += dxx;

                                iyc1Top6++;
                            }

                        }

                        if (arrIntensity[i][iCharge + 1] > 10 * Double.MIN_VALUE) {
                            sumXYMatchPredict += dxy;
                            sumYYMatchPredict += dyy;
                            sumXXMatchPredict += dxx;


                        }

                        //debug
                        if (msOneTrail.getId()==157058 && pep.composition.equals("VVAGVAAALAHK") && ms2Int!=null
                                && ms2Int.arrallTop6IntensityPos!=null){
                            if(Utils.iDebugoutput > 0)
                            {
                                try {
                                    Utils.bufferedWriterDebugInfo.write("rt2:"+ yIonMatchTrail!=null && yIonMatchTrail.length>0 && yIonMatchTrail[currentPosWithMatch]!=null? ""+ yIonMatchTrail[currentPosWithMatch].getRtApex() :"0"+"\n");
                                    Utils.bufferedWriterDebugInfo.write("yi:"+ i+"\n");
                                    Utils.bufferedWriterDebugInfo.write("yIonMatchTrail[currentPosWithMatch].getPeakArea() :"+yIonMatchTrail[currentPosWithMatch].getPeakArea() +"\n");
                                    Utils.bufferedWriterDebugInfo.write("arrIntensity[i][iCharge + 1]:"+arrIntensity[i][iCharge + 1]+"\n");


                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }


                            }

                        }

                    }
//                    else
//                    {
//                        if(debugOutput++<10) {
//                            System.out.println(" no arrIntensity "+pep.composition+"  "+msOneTrail.getZ());
//                        }
//                    }
                } else {
                    if (arrIntensity != null) {
                        double dyy = arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
//                        System.out.println("   ynotin----"+i+"  "+pep.composition+"   "+ msOneTrail.getZ());

//                        System.out.println("arrIntensity:" + JSONWriter.valueToString(arrIntensity));
//                        sumXY += Utils.NoneMatchedIntensity * arrIntensity[i][iCharge + 1];
//                        sumXX += Utils.NoneMatchedIntensity * Utils.NoneMatchedIntensity;

                        sumYY += dyy;
                        dysumYY += dyy;

                        if(ms2Int!=null && ms2Int.arrallTop6IntensityPos!=null) {
                            if (itopall < Utils.iNumPredictIntensityTopN && ms2Int.arrallTop6IntensityCount[iCharge + 1] > 0) {
                                if (ms2Int.arrallTop6IntensityPos[itopall] ==  i ) {
                                    dallTop6sumYY += dyy;
                                    itopall++;

                                }
                            }


                            if (iCharge == 1 && iallc1Top6 < Utils.iNumPredictIntensityTopN && ms2Int.arrallTop6IntensityCharge1 <= iallc1Top6) {
                                if (ms2Int.arrallTop6IntensityCharge1Pos[iallc1Top6] ==  i ) {
                                    dcharget1Top6sumYY += dyy;
                                    iallc1Top6++;
                                }
                            }


                            if (iCharge == 1 && iyc1Top6 < Utils.iNumPredictIntensityTopN && ms2Int.arryTop6IntensityPosCharge1[iyc1Top6] ==  i ) {
                                dyTop6sumYY += dyy;
                                iyc1Top6++;
                            }
                        }

                        //debug
                        if (msOneTrail.getId()==157058 && pep.composition.equals("VVAGVAAALAHK") && ms2Int!=null
                                && ms2Int.arrallTop6IntensityPos!=null){
                            if(Utils.iDebugoutput > 0)
                            {
                                try {
                                    Utils.bufferedWriterDebugInfo.write("yi:"+ i+"\n");
                                    Utils.bufferedWriterDebugInfo.write("arrIntensity[i][iCharge + 1]:"+arrIntensity[i][iCharge + 1]+"\n");


                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }


                            }

                        }

                    }

                }

            }
            //建立特征--不考虑跨rt
            if(ibCountFragmentTypeMatched == 6)//1+2+3 H2O+NH3+none
            {
                ibAllMatchedWithloss ++;
            }else if(ibCountFragmentTypeMatched == 5)//2+3 NH3+none
            {
                ibMatchNormalAndNH3loss++;
            }else if(ibCountFragmentTypeMatched == 4)//1+3 H2O+none
            {
                ibMatchNormalAndH2Oloss ++;
            }
            if(iyCountFragmentTypeMatched == 6)//1+2+3 H2O+NH3+none
            {
                iyAllMatchedWithloss ++;
            }else if(iyCountFragmentTypeMatched == 5)//2+3 NH3+none
            {
                iyMatchNormalAndNH3loss++;
            }else if(iyCountFragmentTypeMatched == 4)//1+3 H2O+none
            {
                iyMatchNormalAndH2Oloss ++;
            }
        }
//            if(msOneTrail.getId()==24474)
//                System.out.println("ok");
        bSum[iCharge-1] = bipos.size();
        ySum[iCharge-1] = yipos.size();
        bySum = bySum + bSum[iCharge-1]+ySum[iCharge-1];

        if (FastMath.abs(dbsumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dbsumYY) < 10 * Double.MIN_VALUE) {
            dBionPearsonCorrelationWithPredictMSMS[iCharge-1] =  0.0;
        }else
        {
            dBionPearsonCorrelationWithPredictMSMS[iCharge-1] = dbsumXY/FastMath.sqrt(dbsumXX*dbsumYY);
        }

        if (FastMath.abs(dysumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dysumYY) < 10 * Double.MIN_VALUE) {
            dYionPearsonCorrelationWithPredictMSMS[iCharge-1] =  0.0;
        }else
        {
            dYionPearsonCorrelationWithPredictMSMS[iCharge-1] = dysumXY/FastMath.sqrt(dysumXX*dysumYY);
        }

        if (bySum > Utils.thresholdFragmentCount/*0*/) {
            if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance) {
                dMatch += dMatch * log10(bySum) / log10(lMSTwotrails.size());
//                dMatch = dMatch * bySum / lMSTwotrails.size();

            }
//                strb += "|B"+iCharge+":" + bipos.stream().map(x -> x + 1).collect(Collectors.toList()).toString();
//                stry += "|Y"+iCharge+":" + yipos.stream().map(x -> len - x).collect(Collectors.toList()).toString();

//            dMatch = dMatch * log10(bySum) / log10(2*pep.b_ions.length);
//            dMatch = dMatch * log10(bySum);


        } else {
            dMatch = 0.0;//未匹配到
        }
    }

    if (bySum >  Utils.thresholdFragmentCount/*0*/) {
        String strDecoy = "target";
//                if(pep.id.contains("Decoy"))//DeBruijnsp
        if (pep.id.contains("DeBruijn"))//>DeBruijnsp
        {
            strDecoy = "decoy";
        }
        String mstwortapex = "" + lMSTwotrails.get(0).getRtApex();
        matchPep = new MatchPeptideWithMultiCharge();
        matchPep.setValues(
                dMatch, pep.id, strDecoy, pep.composition, pep.mass, msOneTrail,
                mstwortapex,  bSum, ySum,  dMatch, strb
                , stry, bipos, yipos);
        matchPep.setPeakAreaMSErrorPeakHalfRankWithMultiCharge(dlogBionPeakAreaSum, dlogYionPeakAreaSum,
                dBionMassErrorSum, dYionMassErrorSum, dBionPeakHalfSum, dYionPeakHalfSum,
                dBionRetentionTimeErrorSum, dYionRetentionTimeErrorSum, pep.dMutationRate);

        matchPep.setConsetivePeakWithMultiCharge(iBConsective, iYConsective);
        matchPep.setPep(pep);

        matchPep.setMapBionMS2Trail(mapBionTrail);
        matchPep.setMapYionMS2Trail(mapYionTrail);
//            matchPep.mapCurRTBionMS2Trail = mapCurRTBionMS2Trail;
//            matchPep.mapCurRTYionMS2Trail = mapCurRTYionMS2Trail;

        matchPep.arr_dBionPearsonCorrelationWithPredictMSMS = dBionPearsonCorrelationWithPredictMSMS;
        matchPep.arr_dYionPearsonCorrelationWithPredictMSMS = dYionPearsonCorrelationWithPredictMSMS;

        //GET THE MAXINUM LOSS COUNT
        if(arrbloss[0]+arrbloss[1]+arrbloss[2]+arryloss[0]+arryloss[1]+arryloss[2]<
                ibAllMatchedWithloss+ibMatchNormalAndNH3loss+ibMatchNormalAndH2Oloss+iyAllMatchedWithloss+iyMatchNormalAndNH3loss+iyMatchNormalAndH2Oloss)
        {
            arrbloss[0] = ibAllMatchedWithloss;
            arrbloss[1] = ibMatchNormalAndNH3loss;
            arrbloss[2] = ibMatchNormalAndH2Oloss;
            arryloss[0] = iyAllMatchedWithloss;
            arryloss[1] = iyMatchNormalAndNH3loss;
            arryloss[2] = iyMatchNormalAndH2Oloss;
        }

        matchPep.ibAllMatchedWithloss    = arrbloss[0];
        matchPep.ibMatchNormalAndNH3loss = arrbloss[1];
        matchPep.ibMatchNormalAndH2Oloss = arrbloss[2];
        matchPep.iyAllMatchedWithloss    = arryloss[0];
        matchPep.iyMatchNormalAndNH3loss = arryloss[1];
        matchPep.iyMatchNormalAndH2Oloss = arryloss[2];

        //B 1 ION 60%
        if (FastMath.abs(db60sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(db60sumYY) < 10 * Double.MIN_VALUE) {
            matchPep.dBion60PearsonCorrelationWithPredictMSMS= 0.0;
        }else
        {
            matchPep.dBion60PearsonCorrelationWithPredictMSMS= db60sumXY/FastMath.sqrt(db60sumXX*db60sumYY);


        }



        //matched and predict cos
        if (FastMath.abs(sumXXMatchPredict) < 10 * Double.MIN_VALUE || FastMath.abs(sumYYMatchPredict) < 10 * Double.MIN_VALUE) {
            matchPep.dCosWithPredictMSMSMatchedPredict= 0.0;
        }else
        {
            matchPep.dCosWithPredictMSMSMatchedPredict=  sumXYMatchPredict/FastMath.sqrt(sumXXMatchPredict*sumYYMatchPredict);
        }

        if (FastMath.abs(sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(sumYY) < 10 * Double.MIN_VALUE) {
            matchPep.dPearsonCorrelationWithPredictMSMS= 0.0;
        }else
        {
            matchPep.dPearsonCorrelationWithPredictMSMS= sumXY/FastMath.sqrt(sumXX*sumYY);
//                if(debugOutput++<100) {
//                    System.out.println("In matchPepCos:" + matchPep.dPearsonCorrelationWithPredictMSMS);
//                }

        }

        if (FastMath.abs(dallTop6sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dallTop6sumYY) < 10 * Double.MIN_VALUE) {
            matchPep.dCosWithPredictMSMSMatchedTop6AllPredict= 0.0;
        }else
        {
            matchPep.dCosWithPredictMSMSMatchedTop6AllPredict=  dallTop6sumXY/FastMath.sqrt(dallTop6sumXX*dallTop6sumYY);
        }


        if (FastMath.abs(dcharget1Top6sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dcharget1Top6sumYY) < 10 * Double.MIN_VALUE) {
            matchPep.dCosWithPredictMSMSMatchedTop6Charge1Predict= 0.0;
        }else
        {
            matchPep.dCosWithPredictMSMSMatchedTop6Charge1Predict=  dcharget1Top6sumXY/FastMath.sqrt(dcharget1Top6sumXX*dcharget1Top6sumYY);
        }

        if (FastMath.abs(dbTop6sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dbTop6sumYY) < 10 * Double.MIN_VALUE) {
            matchPep.dCosWithPredictMSMSMatchedTop6BIonCharge1Predict= 0.0;
        }else
        {
            matchPep.dCosWithPredictMSMSMatchedTop6BIonCharge1Predict=  dbTop6sumXY/FastMath.sqrt(dbTop6sumXX*dbTop6sumYY);
        }

        if (FastMath.abs(dyTop6sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dyTop6sumYY) < 10 * Double.MIN_VALUE) {
            matchPep.dCosWithPredictMSMSMatchedTop6YIonCharge1Predict= 0.0;
        }else
        {
            matchPep.dCosWithPredictMSMSMatchedTop6YIonCharge1Predict=  dyTop6sumXY/FastMath.sqrt(dyTop6sumXX*dyTop6sumYY);
        }


//debug
        if (msOneTrail.getId()==157058 && pep.composition.equals("VVAGVAAALAHK") && ms2Int!=null
                && ms2Int.arrallTop6IntensityPos!=null){
            if(Utils.iDebugoutput > 0)
            {
                try {
                    Utils.bufferedWriterDebugInfo.write("bion:"+ JSONWriter.valueToString(pep.b_ions)+"\n");
                    Utils.bufferedWriterDebugInfo.write("yion:"+ JSONWriter.valueToString(pep.y_ions)+"\n");


                } catch (IOException e) {
                    throw new RuntimeException(e);
                }


            }

        }


    }
    //       return dMatch;
    return  matchPep ;

}
    public  static  MatchPeptideWithMultiCharge PepMatchIonCountCrossRTMultiChargeFragmentLoss(Peptide pep, List<MSTwoTrail> lMSTwotrails, double maxPeakArea,
                                                          Enums.ScoreMode scoreMode, BufferedWriter writer, int icount,
                                                          MSOneTrail msOneTrail/*double ms1rt, double ms1Mass, long ms1id,double ms1qualityScore*/,
                                                          MSTwoTrail[] bIonMatchTrail,MSTwoTrail[] yIonMatchTrail,int[] arrbloss,int[] arryloss,
                                                                                   double[][] arrIntensity) {
        List<Integer> bpos = null;
        List<Integer> ypos = null;
        List<Integer> bNormalpos  = Collections.synchronizedList(new ArrayList<>());
        List<Integer> yNormalpos  = Collections.synchronizedList(new ArrayList<>());
        List<Integer> bH2OLosspos = Collections.synchronizedList(new ArrayList<>());
        List<Integer> yH2OLosspos = Collections.synchronizedList(new ArrayList<>());
        List<Integer> bNH3Losspos = Collections.synchronizedList(new ArrayList<>());
        List<Integer> yNH3Losspos = Collections.synchronizedList(new ArrayList<>());

        List<Integer> bipos = Collections.synchronizedList(new ArrayList<>());
        List<Integer> yipos = Collections.synchronizedList(new ArrayList<>());

//        List<Double> bMatchDistance =  Collections.synchronizedList(new ArrayList<>());
//        List<Double> yMatchDistance = Collections.synchronizedList(new ArrayList<>());

        double dMatch = 1.0;
        double[] dlogBionPeakAreaSum = new double[Utils.iMS2Charge];
        double[] dlogYionPeakAreaSum = new double[Utils.iMS2Charge];
        double[] dBionMassErrorSum = new double[Utils.iMS2Charge];
        double[] dYionMassErrorSum = new double[Utils.iMS2Charge];
        double[] dBionRetentionTimeErrorSum = new double[Utils.iMS2Charge];
        double[] dYionRetentionTimeErrorSum = new double[Utils.iMS2Charge];
        double[] dBionPeakHalfSum = new double[Utils.iMS2Charge];
        double[] dYionPeakHalfSum = new double[Utils.iMS2Charge];
        int[] iBConsective = new int[Utils.iMS2Charge];
        int[] iYConsective = new int[Utils.iMS2Charge];
        String strb = "";
        String stry = "";
        int[] bSum = new int[Utils.iMS2Charge];
        int[] ySum = new int[Utils.iMS2Charge];
        int bySum = 0;
        int iBYComplementaryTrailCount = 0;

        double sumXX = 0d;
        double sumYY = 0d;
        double sumXY = 0d;

        double sumXXMatchPredict = 0d;
        double sumYYMatchPredict = 0d;
        double sumXYMatchPredict = 0d;

        double db60sumXX = 0d;
        double db60sumYY = 0d;
        double db60sumXY = 0d;

        TreeMap<Integer, LinkedList<MSTwoTrail>>  mapBionTrail = null;
        TreeMap<Integer, LinkedList<MSTwoTrail>>  mapYionTrail = null;

//        TreeMap<Integer, MSTwoTrail>  mapCurRTBionMS2Trail = null;
//        TreeMap<Integer, MSTwoTrail>  mapCurRTYionMS2Trail = null;

        double[] dBionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];
        double[] dYionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];

        int ibAllMatchedWithloss = 0;
        int ibMatchNormalAndH2Oloss = 0;
        int ibMatchNormalAndNH3loss = 0;
        int iyAllMatchedWithloss = 0;
        int iyMatchNormalAndH2Oloss = 0;
        int iyMatchNormalAndNH3loss = 0;

        if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances || Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
        {
            dMatch = 0.0;//add all abundance
        }
        int len = pep.b_ions.length;
        MatchPeptideWithMultiCharge matchPep = null;

//        int iCurrentPos = 0;
        for(int iCharge=1;iCharge<=Utils.iMS2Charge;iCharge++) {
//            bpos.clear();
//            ypos.clear();
            bNormalpos.clear();
            yNormalpos.clear();
            bH2OLosspos.clear();
            yH2OLosspos.clear();
            bNH3Losspos.clear();
            yNH3Losspos.clear();

            bipos.clear();
            yipos.clear();


            int ibposFront = 0;
            int ibposBack = 0;
            int iyposFront = 0;
            int iyposBack = 0;

            int ibnearpos = 0;
            int iynearpos = 0;

            int bMatchSize;
            int yMatchSize;

            double dbsumXX = 0d;
            double dbsumYY = 0d;
            double dbsumXY = 0d;

            double dysumXX = 0d;
            double dysumYY = 0d;
            double dysumXY = 0d;




            for (int i = 0; i < len -1; ++i) {// b y  ion less than len of peptides

                //发现fragment loos有两种策略，一是将丢失的离子看成是一个正常的可匹配的离子峰，二是将其作为已匹配的正常离子峰的附加属性（若匹配的正常的离子再检查他有没有附加的丢失水和胺的离子峰并统计其个数）
                //下面的代码是策略一
                int ibCountFragmentTypeMatched = 0;
                int iyCountFragmentTypeMatched = 0;

                for (int iFragmentType = 0; iFragmentType < Utils.iMS2FragmentTypeCount; iFragmentType++) {

                    double b_ionsMZ = (pep.b_ions[i] + AminoAcid.ProtonMass * (iCharge - 1)) / iCharge;//由于b离子原来加过一个离子，对后面的离子只需要加入一个质子质量 即可 mass to mz 该公式化简是 mass/z + AminoAcid.ProtonMass;
                    double y_ionsMZ = (pep.y_ions[len - i - 1] + AminoAcid.ProtonMass * (iCharge - 1)) / iCharge;//由于y离子原来加过一个离子，对后面的离子只需要加入一个质子质量 即可
                    if(Utils.iMS2FragmentTypeCount==3) {//考虑loss
                        //从丢失的 大小开始，越大越在前面，mz越小，这样search时候先搜小的然后搜大的
                        if (iFragmentType == 0) {
                            bpos = bH2OLosspos;
                            ypos = yH2OLosspos;
                            b_ionsMZ = b_ionsMZ - (Utils.H2OMass) / iCharge;
                            y_ionsMZ = y_ionsMZ - (Utils.H2OMass) / iCharge;
                        } else if (iFragmentType == 1) {
                            bpos = bNH3Losspos;
                            ypos = yNH3Losspos;
                            b_ionsMZ = b_ionsMZ - (Utils.NH3Mass) / iCharge;
                            y_ionsMZ = y_ionsMZ - (Utils.NH3Mass) / iCharge;
                        } else {
                            bpos = bNormalpos;
                            ypos = yNormalpos;
                        }
                    }else if(Utils.iMS2FragmentTypeCount==1)//不考虑loss
                    {
                        bpos = bNormalpos;
                        ypos = yNormalpos;
                    }
                    bMatchSize = bpos.size();
                    yMatchSize = ypos.size();

//                double minBIonDistance = Utils.thresholdPPM;
//                double minYIonDistance = Utils.thresholdPPM;

                    double[][] arr_lftIsotope = new double[Utils.iMS2Charge][2];

                    double[][] arr_rightIsotope = new double[Utils.iMS2Charge][2];
                    for (int ic = 0; ic < Utils.iMS2Charge; ic++) {
                        arr_lftIsotope[ic][0] = b_ionsMZ * (1 - Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                        arr_lftIsotope[ic][1] = b_ionsMZ * (1 + Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                        arr_rightIsotope[ic][0] = b_ionsMZ * (1 - Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                        arr_rightIsotope[ic][1] = b_ionsMZ * (1 + Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                    }

                    double bIonMinusMass = b_ionsMZ * (1 - Utils.thresholdPPM);
                    double bIonMaxMass = b_ionsMZ * (1 + Utils.thresholdPPM);

                    int iBionCurCount = 0;
                    int iYionCurCount = 0;


                    //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
                    if ((bIonMaxMass > (lMSTwotrails.get(0).getMzApex()))
                            && (bIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size() - 1).getMzApex()))) {
                        ibposFront = binarySearch0(lMSTwotrails, ibposFront, lMSTwotrails.size(), bIonMinusMass);//全新查找
                        ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposFront);

                        ///更新bpos里面的内容
//                    iBionCurCount = getiIonCurCount(lMSTwotrails, ibnearpos, bIonMaxMass, bpos);
//                    iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass,b_ionsPreIsotopeMaxMZ,b_ionsPreIsotopeMinMZ
//                            ,arr_rightIsotope,iCharge, bpos);
                        iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass, arr_lftIsotope
                                , arr_rightIsotope, iCharge, bpos);

                        if (bpos.size() > 0) {
                            ibposFront = bpos.get(bpos.size() - 1);
                        }
                        if (ibposFront < 0) {
                            ibposFront = ibnearpos;
                        }
                    }


//                double y_ionsPreIsotopeMinMZ = y_ionsMZ * (1 -  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / iCharge;
//                double y_ionsPreIsotopeMaxMZ = y_ionsMZ * (1 +  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / iCharge;


                    for (int ic = 0; ic < Utils.iMS2Charge; ic++) {
                        arr_lftIsotope[ic][0] = y_ionsMZ * (1 - Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                        arr_lftIsotope[ic][1] = y_ionsMZ * (1 + Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                        arr_rightIsotope[ic][0] = y_ionsMZ * (1 - Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);
                        arr_rightIsotope[ic][1] = y_ionsMZ * (1 + Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic + 1);

                    }


                    double yIonMinusMass = y_ionsMZ * (1 - Utils.thresholdPPM);
                    double yIonMaxMass = y_ionsMZ * (1 + Utils.thresholdPPM);


                    //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
                    if ((yIonMaxMass > (lMSTwotrails.get(0).getMzApex()))
                            && (yIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size() - 1).getMzApex()))) {
                        //根据前面比较的结果，search 对应的y离子

                        iyposFront = binarySearch0(lMSTwotrails, iyposFront, lMSTwotrails.size(), yIonMinusMass);
                        iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);

//                    iYionCurCount = getiIonCurCount(lMSTwotrails, iynearpos, yIonMaxMass, ypos);
//                    iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass,y_ionsPreIsotopeMaxMZ,y_ionsPreIsotopeMinMZ
//                            ,arr_rightIsotope,iCharge, ypos);
                        iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass, arr_lftIsotope
                                , arr_rightIsotope, iCharge, ypos);


                        if (ypos.size() > 0) {
                            iyposFront = ypos.get(ypos.size() - 1);
                        }
                        if (iyposFront < 0) {
                            iyposFront = iynearpos;
                        }
                    }


                    //判断单个bion with ms2 trail匹配是否成功
                    //判断连续ion
                    //大于某阈值后的连续离子，大于ms2的中平均数？中位数？大于peptide的质量的一半？
                    //根据ms2 trail数据的散列度设置阈值？

                    //用于记录哪个离子峰匹配，包含：loss类型、电荷、离子峰编号
                    //framentType+icharge+iNumber,没有loss排前面
                    int currentPosWithMatch = (Utils.iMS2FragmentTypeCount-iFragmentType-1) * (Utils.iMS2Charge * len) +(iCharge - 1) * len + i  ;
                    if (bpos.size() > bMatchSize) {
                        ibCountFragmentTypeMatched += (iFragmentType + 1);//建立特征--不考虑跨rt
//                        if (mapCurRTBionMS2Trail == null) {
//                            mapCurRTBionMS2Trail = new TreeMap<>();
//                        }
//                        mapCurRTBionMS2Trail.put(currentPosWithMatch + 1, lMSTwotrails.get(ibposFront));//其他地方没有用到
                        if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                        {
                            bIonMatchTrail[currentPosWithMatch] = lMSTwotrails.get(ibposFront);
                        } else//匹配多个MS2RT
                        {
                            bIonMatchTrail[currentPosWithMatch] = CompateTwoTrail(lMSTwotrails.get(ibposFront), bIonMatchTrail[currentPosWithMatch], b_ionsMZ);
                        }
                    } else {
                        if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                        {
                            bIonMatchTrail[currentPosWithMatch] = null;
                        }
                    }

                    if (bIonMatchTrail[currentPosWithMatch] != null)//bpos.size()>bMatchSize)
                    {
                        if ((bipos.size() > 0) && (bipos.get(bipos.size() - 1) == currentPosWithMatch - 1)) {
                            iBConsective[iCharge - 1]++;
                        }

                        bipos.add(currentPosWithMatch);

                        if (mapBionTrail == null) {
                            mapBionTrail = new TreeMap<>();
                        }
                        LinkedList<MSTwoTrail> listMS2Trail = mapBionTrail.get(currentPosWithMatch);
                        if (listMS2Trail == null) {
                            listMS2Trail = new LinkedList<MSTwoTrail>();
                        }
                        listMS2Trail.add(bIonMatchTrail[currentPosWithMatch]);
                        mapBionTrail.put(currentPosWithMatch, listMS2Trail);

                        dMatch = getdMatchIonScore(bIonMatchTrail[currentPosWithMatch], bIonMaxMass, dMatch, maxPeakArea, iBionCurCount);

                        dMatch = getdMatchMultiIonScore(lMSTwotrails, bipos, bpos, dMatch);

                        //用于输出中间值

                        dlogBionPeakAreaSum[iCharge - 1] += Math.log10(bIonMatchTrail[currentPosWithMatch].getPeakArea());
                        dBionMassErrorSum[iCharge - 1] += 1.0 - Math.pow((((bIonMatchTrail[currentPosWithMatch].getMzApex() - b_ionsMZ) / (bIonMatchTrail[currentPosWithMatch].getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                        dBionPeakHalfSum[iCharge - 1] += Math.max(Math.log10(50.0 / bIonMatchTrail[currentPosWithMatch].getPeakAreaHalfRank()), 0.0);
                        dBionRetentionTimeErrorSum[iCharge - 1] += 1 - Math.pow((bIonMatchTrail[currentPosWithMatch].getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
                        if (arrIntensity != null) {
//                        System.out.println("B arrIntensity "+pep.composition+"  "+msOneTrail.getZ());

                            sumXY += bIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge - 1];
                            sumYY += arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];
                            sumXX += bIonMatchTrail[currentPosWithMatch].getPeakArea() * bIonMatchTrail[currentPosWithMatch].getPeakArea();

                            dbsumXY += bIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge - 1];
                            dbsumYY += arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];
                            dbsumXX += bIonMatchTrail[currentPosWithMatch].getPeakArea() * bIonMatchTrail[currentPosWithMatch].getPeakArea();

                            if (iCharge == 1 & i < len * 0.6) {
                                db60sumXY += bIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge - 1];
                                db60sumYY += arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];
                                db60sumXX += bIonMatchTrail[currentPosWithMatch].getPeakArea() * bIonMatchTrail[currentPosWithMatch].getPeakArea();

                            }

                            if (arrIntensity[i][iCharge - 1] > 10 * Double.MIN_VALUE) {
                                sumXYMatchPredict += bIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge - 1];
                                sumYYMatchPredict += arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];
                                sumXXMatchPredict += bIonMatchTrail[currentPosWithMatch].getPeakArea() * bIonMatchTrail[currentPosWithMatch].getPeakArea();

                            }
                        }
//                    else
//                    {
//                        if(debugOutput++<10) {
//                            System.out.println(" no arrIntensity "+pep.composition+"  "+msOneTrail.getZ());
//                        }
//                    }

                    } else {

                        if (arrIntensity != null) {
//                        System.out.println("   bnotin----"+i+"  "+pep.composition+"   "+ msOneTrail.getZ());

//                        System.out.println("arrIntensity:" + JSONWriter.valueToString(arrIntensity));
//                        sumXY += Utils.NoneMatchedIntensity * arrIntensity[i][iCharge -1 ];
//                        sumXX += Utils.NoneMatchedIntensity * Utils.NoneMatchedIntensity;
                            sumYY += arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];//做cosine时，没有匹配Y值（不是y离子）也要加

                            dbsumYY += arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];//做cosine时，没有匹配Y值（不是y离子）也要加

                            if (iCharge == 1 & i < len * 0.6) {
                                db60sumYY += arrIntensity[i][iCharge - 1] * arrIntensity[i][iCharge - 1];//做cosine时，没有匹配Y值（不是y离子）也要加

                            }
                        }
                    }

                    //用于记录哪个离子峰匹配，包含：loss类型、电荷、离子峰编号
                    //framentType+icharge+iNumber,没有loss排前面
                     currentPosWithMatch = (Utils.iMS2FragmentTypeCount-iFragmentType-1) * (Utils.iMS2Charge * len) + (iCharge - 1) * len + len - i - 1  ;

                    if (ypos.size() > yMatchSize) {
                        iyCountFragmentTypeMatched += (iFragmentType + 1);//建立特征--不考虑跨rt

//                        if (mapCurRTYionMS2Trail == null) {
//                            mapCurRTYionMS2Trail = new TreeMap<>();
//                        }
//                        mapCurRTYionMS2Trail.put((iCharge - 1) * len + i + 1, lMSTwotrails.get(iyposFront));//其他地方没有用到
                        if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                        {
                            yIonMatchTrail[currentPosWithMatch] = lMSTwotrails.get(iyposFront);
                        } else//匹配多个MS2RT
                        {
                            yIonMatchTrail[currentPosWithMatch] = CompateTwoTrail(lMSTwotrails.get(iyposFront), yIonMatchTrail[currentPosWithMatch], y_ionsMZ);
                        }
                    } else {
                        if (Utils.bzMatchOneRT)//只匹配一个ms2RT
                        {
                            yIonMatchTrail[currentPosWithMatch] = null;
                        }
                    }

                    if (yIonMatchTrail[currentPosWithMatch] != null) {
                        if ((yipos.size() > 0) && (yipos.get(yipos.size() - 1) == (iCharge - 1) * len + len - i)) {
                            iYConsective[iCharge - 1]++;
                        }


                        yipos.add(currentPosWithMatch);

                        if (mapYionTrail == null) {
                            mapYionTrail = new TreeMap<>();
                        }
                        LinkedList<MSTwoTrail> listMS2Trail = mapYionTrail.get(currentPosWithMatch);
                        if (listMS2Trail == null) {
                            listMS2Trail = new LinkedList<MSTwoTrail>();
                        }
                        listMS2Trail.add(yIonMatchTrail[currentPosWithMatch]);
                        mapYionTrail.put(currentPosWithMatch, listMS2Trail);


                        dMatch = getdMatchIonScore(yIonMatchTrail[currentPosWithMatch], yIonMaxMass, dMatch, maxPeakArea, iYionCurCount);
                        dMatch = getdMatchMultiIonScore(lMSTwotrails, yipos, ypos, dMatch);

                        dlogYionPeakAreaSum[iCharge - 1] += Math.log10(yIonMatchTrail[currentPosWithMatch].getPeakArea());
                        dYionMassErrorSum[iCharge - 1] += 1.0 - Math.pow((((yIonMatchTrail[currentPosWithMatch].getMzApex() - y_ionsMZ) / (yIonMatchTrail[currentPosWithMatch].getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                        dYionPeakHalfSum[iCharge - 1] += Math.max(Math.log10(50.0 / yIonMatchTrail[currentPosWithMatch].getPeakAreaHalfRank()), 0.0);
                        dYionRetentionTimeErrorSum[iCharge - 1] += 1 - Math.pow((yIonMatchTrail[currentPosWithMatch].getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
                        if (arrIntensity != null) {
//                        System.out.println("Y arrIntensity "+pep.composition+"  "+msOneTrail.getZ());

                            sumXY += yIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge + 1];
                            sumYY += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                            sumXX += yIonMatchTrail[currentPosWithMatch].getPeakArea() * yIonMatchTrail[currentPosWithMatch].getPeakArea();

                            dysumXY += yIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge + 1];
                            dysumYY += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                            dysumXX += yIonMatchTrail[currentPosWithMatch].getPeakArea() * yIonMatchTrail[currentPosWithMatch].getPeakArea();

                            if (arrIntensity[i][iCharge + 1] > 10 * Double.MIN_VALUE) {
                                sumXYMatchPredict += yIonMatchTrail[currentPosWithMatch].getPeakArea() * arrIntensity[i][iCharge + 1];
                                sumYYMatchPredict += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                                sumXXMatchPredict += yIonMatchTrail[currentPosWithMatch].getPeakArea() * yIonMatchTrail[currentPosWithMatch].getPeakArea();


                            }

                        }
//                    else
//                    {
//                        if(debugOutput++<10) {
//                            System.out.println(" no arrIntensity "+pep.composition+"  "+msOneTrail.getZ());
//                        }
//                    }
                    } else {
                        if (arrIntensity != null) {
//                        System.out.println("   ynotin----"+i+"  "+pep.composition+"   "+ msOneTrail.getZ());

//                        System.out.println("arrIntensity:" + JSONWriter.valueToString(arrIntensity));
//                        sumXY += Utils.NoneMatchedIntensity * arrIntensity[i][iCharge + 1];
//                        sumXX += Utils.NoneMatchedIntensity * Utils.NoneMatchedIntensity;

                            sumYY += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                            dysumYY += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                        }

                    }

                }
                //建立特征--不考虑跨rt
                if(ibCountFragmentTypeMatched == 6)//1+2+3 H2O+NH3+none
                {
                    ibAllMatchedWithloss ++;
                }else if(ibCountFragmentTypeMatched == 5)//2+3 NH3+none
                {
                    ibMatchNormalAndNH3loss++;
                }else if(ibCountFragmentTypeMatched == 4)//1+3 H2O+none
                {
                    ibMatchNormalAndH2Oloss ++;
                }
                if(iyCountFragmentTypeMatched == 6)//1+2+3 H2O+NH3+none
                {
                    iyAllMatchedWithloss ++;
                }else if(iyCountFragmentTypeMatched == 5)//2+3 NH3+none
                {
                    iyMatchNormalAndNH3loss++;
                }else if(iyCountFragmentTypeMatched == 4)//1+3 H2O+none
                {
                    iyMatchNormalAndH2Oloss ++;
                }
            }
//            if(msOneTrail.getId()==24474)
//                System.out.println("ok");
            bSum[iCharge-1] = bipos.size();
            ySum[iCharge-1] = yipos.size();
            bySum = bySum + bSum[iCharge-1]+ySum[iCharge-1];

            if (FastMath.abs(dbsumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dbsumYY) < 10 * Double.MIN_VALUE) {
                dBionPearsonCorrelationWithPredictMSMS[iCharge-1] =  0.0;
            }else
            {
                dBionPearsonCorrelationWithPredictMSMS[iCharge-1] = dbsumXY/FastMath.sqrt(dbsumXX*dbsumYY);
            }

            if (FastMath.abs(dysumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dysumYY) < 10 * Double.MIN_VALUE) {
                dYionPearsonCorrelationWithPredictMSMS[iCharge-1] =  0.0;
            }else
            {
                dYionPearsonCorrelationWithPredictMSMS[iCharge-1] = dysumXY/FastMath.sqrt(dysumXX*dysumYY);
            }

            if (bySum > Utils.thresholdFragmentCount/*0*/) {
                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance) {
                    dMatch += dMatch * log10(bySum) / log10(lMSTwotrails.size());
//                dMatch = dMatch * bySum / lMSTwotrails.size();

                }
//                strb += "|B"+iCharge+":" + bipos.stream().map(x -> x + 1).collect(Collectors.toList()).toString();
//                stry += "|Y"+iCharge+":" + yipos.stream().map(x -> len - x).collect(Collectors.toList()).toString();

//            dMatch = dMatch * log10(bySum) / log10(2*pep.b_ions.length);
//            dMatch = dMatch * log10(bySum);


            } else {
                dMatch = 0.0;//未匹配到
            }
        }

        if (bySum >  Utils.thresholdFragmentCount/*0*/) {
            String strDecoy = "target";
//                if(pep.id.contains("Decoy"))//DeBruijnsp
            if (pep.id.contains("DeBruijn"))//>DeBruijnsp
            {
                strDecoy = "decoy";
            }
            String mstwortapex = "" + lMSTwotrails.get(0).getRtApex();
            matchPep = new MatchPeptideWithMultiCharge();
            matchPep.setValues(
                     dMatch, pep.id, strDecoy, pep.composition, pep.mass, msOneTrail,
                    mstwortapex,  bSum, ySum,  dMatch, strb
                    , stry, bipos, yipos);
            matchPep.setPeakAreaMSErrorPeakHalfRankWithMultiCharge(dlogBionPeakAreaSum, dlogYionPeakAreaSum,
                    dBionMassErrorSum, dYionMassErrorSum, dBionPeakHalfSum, dYionPeakHalfSum,
                    dBionRetentionTimeErrorSum, dYionRetentionTimeErrorSum, pep.dMutationRate);

            matchPep.setConsetivePeakWithMultiCharge(iBConsective, iYConsective);
            matchPep.setPep(pep);

            matchPep.setMapBionMS2Trail(mapBionTrail);
            matchPep.setMapYionMS2Trail(mapYionTrail);
//            matchPep.mapCurRTBionMS2Trail = mapCurRTBionMS2Trail;
//            matchPep.mapCurRTYionMS2Trail = mapCurRTYionMS2Trail;

            matchPep.arr_dBionPearsonCorrelationWithPredictMSMS = dBionPearsonCorrelationWithPredictMSMS;
            matchPep.arr_dYionPearsonCorrelationWithPredictMSMS = dYionPearsonCorrelationWithPredictMSMS;

            //GET THE MAXINUM LOSS COUNT
            if(arrbloss[0]+arrbloss[1]+arrbloss[2]+arryloss[0]+arryloss[1]+arryloss[2]<
                    ibAllMatchedWithloss+ibMatchNormalAndNH3loss+ibMatchNormalAndH2Oloss+iyAllMatchedWithloss+iyMatchNormalAndNH3loss+iyMatchNormalAndH2Oloss)
            {
                arrbloss[0] = ibAllMatchedWithloss;
                arrbloss[1] = ibMatchNormalAndNH3loss;
                arrbloss[2] = ibMatchNormalAndH2Oloss;
                arryloss[0] = iyAllMatchedWithloss;
                arryloss[1] = iyMatchNormalAndNH3loss;
                arryloss[2] = iyMatchNormalAndH2Oloss;
            }

            matchPep.ibAllMatchedWithloss    = arrbloss[0];
            matchPep.ibMatchNormalAndNH3loss = arrbloss[1];
            matchPep.ibMatchNormalAndH2Oloss = arrbloss[2];
            matchPep.iyAllMatchedWithloss    = arryloss[0];
            matchPep.iyMatchNormalAndNH3loss = arryloss[1];
            matchPep.iyMatchNormalAndH2Oloss = arryloss[2];

            //B 1 ION 60%
            if (FastMath.abs(db60sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(db60sumYY) < 10 * Double.MIN_VALUE) {
                matchPep.dBion60PearsonCorrelationWithPredictMSMS= 0.0;
            }else
            {
                matchPep.dBion60PearsonCorrelationWithPredictMSMS= db60sumXY/FastMath.sqrt(db60sumXX*db60sumYY);


            }

            //matched and predict cos
            if (FastMath.abs(sumXXMatchPredict) < 10 * Double.MIN_VALUE || FastMath.abs(sumYYMatchPredict) < 10 * Double.MIN_VALUE) {
                matchPep.dCosWithPredictMSMSMatchedPredict= 0.0;
            }else
            {
                matchPep.dCosWithPredictMSMSMatchedPredict=  sumXYMatchPredict/FastMath.sqrt(sumXXMatchPredict*sumYYMatchPredict);


            }

            if (FastMath.abs(sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(sumYY) < 10 * Double.MIN_VALUE) {
                matchPep.dPearsonCorrelationWithPredictMSMS= 0.0;
            }else
            {
                matchPep.dPearsonCorrelationWithPredictMSMS= sumXY/FastMath.sqrt(sumXX*sumYY);
//                if(debugOutput++<100) {
//                    System.out.println("In matchPepCos:" + matchPep.dPearsonCorrelationWithPredictMSMS);
//                }

            }


        }
        //       return dMatch;
        return  matchPep ;

    }
    public  static  MatchPeptideWithMultiCharge PepMatchIonCountCrossRTMultiCharge(Peptide pep, List<MSTwoTrail> lMSTwotrails, double maxPeakArea,
                                                                                   Enums.ScoreMode scoreMode, BufferedWriter writer, int icount,
                                                                                   MSOneTrail msOneTrail/*double ms1rt, double ms1Mass, long ms1id,double ms1qualityScore*/,
                                                                                   MSTwoTrail[] bIonMatchTrail,MSTwoTrail[] yIonMatchTrail,
                                                                                   double[][] arrIntensity) {
        List<Integer> bpos = Collections.synchronizedList(new ArrayList<>());
        List<Integer> ypos = Collections.synchronizedList(new ArrayList<>());

        List<Integer> bipos = Collections.synchronizedList(new ArrayList<>());
        List<Integer> yipos = Collections.synchronizedList(new ArrayList<>());

//        List<Double> bMatchDistance =  Collections.synchronizedList(new ArrayList<>());
//        List<Double> yMatchDistance = Collections.synchronizedList(new ArrayList<>());

        double dMatch = 1.0;
        double[] dlogBionPeakAreaSum = new double[Utils.iMS2Charge];
        double[] dlogYionPeakAreaSum = new double[Utils.iMS2Charge];
        double[] dBionMassErrorSum = new double[Utils.iMS2Charge];
        double[] dYionMassErrorSum = new double[Utils.iMS2Charge];
        double[] dBionRetentionTimeErrorSum = new double[Utils.iMS2Charge];
        double[] dYionRetentionTimeErrorSum = new double[Utils.iMS2Charge];
        double[] dBionPeakHalfSum = new double[Utils.iMS2Charge];
        double[] dYionPeakHalfSum = new double[Utils.iMS2Charge];
        int[] iBConsective = new int[Utils.iMS2Charge];
        int[] iYConsective = new int[Utils.iMS2Charge];
        String strb = "";
        String stry = "";
        int[] bSum = new int[Utils.iMS2Charge];
        int[] ySum = new int[Utils.iMS2Charge];
        int bySum = 0;
        int iBYComplementaryTrailCount = 0;

        double sumXX = 0d;
        double sumYY = 0d;
        double sumXY = 0d;

        double sumXXMatchPredict = 0d;
        double sumYYMatchPredict = 0d;
        double sumXYMatchPredict = 0d;

        double db60sumXX = 0d;
        double db60sumYY = 0d;
        double db60sumXY = 0d;

        TreeMap<Integer, LinkedList<MSTwoTrail>>  mapBionTrail = null;
        TreeMap<Integer, LinkedList<MSTwoTrail>>  mapYionTrail = null;

        TreeMap<Integer, MSTwoTrail>  mapCurRTBionMS2Trail = null;
        TreeMap<Integer, MSTwoTrail>  mapCurRTYionMS2Trail = null;

        double[] dBionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];
        double[] dYionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];


        if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances || Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
        {
            dMatch = 0.0;//add all abundance
        }
        int len = pep.b_ions.length;
        MatchPeptideWithMultiCharge matchPep = null;

//        int iCurrentPos = 0;
        for(int iCharge=1;iCharge<=Utils.iMS2Charge;iCharge++) {
            bpos.clear();
            ypos.clear();

            bipos.clear();
            yipos.clear();


            int ibposFront = 0;
            int ibposBack = 0;
            int iyposFront = 0;
            int iyposBack = 0;

            int ibnearpos = 0;
            int iynearpos = 0;

            int bMatchSize;
            int yMatchSize;

            double dbsumXX = 0d;
            double dbsumYY = 0d;
            double dbsumXY = 0d;

            double dysumXX = 0d;
            double dysumYY = 0d;
            double dysumXY = 0d;




            for (int i = 0; i < len -1; ++i) {// b y  ion less than len of peptides


                bMatchSize = bpos.size();
                yMatchSize = ypos.size();

//                double minBIonDistance = Utils.thresholdPPM;
//                double minYIonDistance = Utils.thresholdPPM;

                double b_ionsMZ = (pep.b_ions[i] + AminoAcid.ProtonMass * (iCharge - 1)) / iCharge;//由于b离子原来加过一个离子，对后面的离子只需要加入一个质子质量 即可 mass to mz 该公式化简是 mass/z + AminoAcid.ProtonMass;

                double b_ionsLossH2OMZ = b_ionsMZ - (Utils.H2OMass) / iCharge;
                double b_ionsLossNH3MZ = b_ionsMZ - (Utils.NH3Mass) / iCharge;
                double[][] arr_lftIsotope = new double[Utils.iMS2Charge][2];

                double[][] arr_rightIsotope = new double[Utils.iMS2Charge][2];
                for(int ic = 0;ic<Utils.iMS2Charge;ic++)
                {
                    arr_lftIsotope[ic][0] = b_ionsMZ * (1 -  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) /  (ic+1);
                    arr_lftIsotope[ic][1] = b_ionsMZ * (1 +  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) /  (ic+1);

                    arr_rightIsotope[ic][0] = b_ionsMZ * (1 -  Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic+1);
                    arr_rightIsotope[ic][1] = b_ionsMZ * (1 +  Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic+1);

                }

                double bIonMinusMass = b_ionsMZ * (1 - Utils.thresholdPPM);
                double bIonMaxMass = b_ionsMZ * (1 + Utils.thresholdPPM);

                int iBionCurCount = 0;
                int iYionCurCount = 0;


                //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
                if ((bIonMaxMass > (lMSTwotrails.get(0).getMzApex()))
                        && (bIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size() - 1).getMzApex()))) {
                    ibposFront = binarySearch0(lMSTwotrails, ibposFront, lMSTwotrails.size(), bIonMinusMass);//全新查找
                    ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposFront);

                    ///更新bpos里面的内容
//                    iBionCurCount = getiIonCurCount(lMSTwotrails, ibnearpos, bIonMaxMass, bpos);
//                    iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass,b_ionsPreIsotopeMaxMZ,b_ionsPreIsotopeMinMZ
//                            ,arr_rightIsotope,iCharge, bpos);
                    iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass,arr_lftIsotope
                            ,arr_rightIsotope,iCharge, bpos);

                    if (bpos.size() > 0) {
                        ibposFront = bpos.get(bpos.size() - 1);
                    }
                    if (ibposFront < 0) {
                        ibposFront = ibnearpos;
                    }
                }

                double y_ionsMZ = (pep.y_ions[len - i - 1] + AminoAcid.ProtonMass * (iCharge - 1)) / iCharge;//由于y离子原来加过一个离子，对后面的离子只需要加入一个质子质量 即可

//                double y_ionsPreIsotopeMinMZ = y_ionsMZ * (1 -  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / iCharge;
//                double y_ionsPreIsotopeMaxMZ = y_ionsMZ * (1 +  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / iCharge;


                for(int ic = 0;ic<Utils.iMS2Charge;ic++)
                {
                    arr_lftIsotope[ic][0] = y_ionsMZ * (1 -  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic+1);
                    arr_lftIsotope[ic][1] = y_ionsMZ * (1 +  Utils.thresholdPPM) - (Utils.C13_MASS - Utils.C12_MASS) / (ic+1);

                    arr_rightIsotope[ic][0] = y_ionsMZ * (1 -  Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic+1);
                    arr_rightIsotope[ic][1] = y_ionsMZ * (1 +  Utils.thresholdPPM) + (Utils.C13_MASS - Utils.C12_MASS) / (ic+1);

                }


                double yIonMinusMass = y_ionsMZ * (1 - Utils.thresholdPPM);
                double yIonMaxMass = y_ionsMZ * (1 + Utils.thresholdPPM);


                //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
                if ((yIonMaxMass > (lMSTwotrails.get(0).getMzApex()))
                        && (yIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size() - 1).getMzApex()))) {
                    //根据前面比较的结果，search 对应的y离子

                    iyposFront = binarySearch0(lMSTwotrails, iyposFront, lMSTwotrails.size(), yIonMinusMass);
                    iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);

//                    iYionCurCount = getiIonCurCount(lMSTwotrails, iynearpos, yIonMaxMass, ypos);
//                    iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass,y_ionsPreIsotopeMaxMZ,y_ionsPreIsotopeMinMZ
//                            ,arr_rightIsotope,iCharge, ypos);
                    iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass,arr_lftIsotope
                            ,arr_rightIsotope,iCharge, ypos);


                    if (ypos.size() > 0) {
                        iyposFront = ypos.get(ypos.size() - 1);
                    }
                    if (iyposFront < 0) {
                        iyposFront = iynearpos;
                    }
                }


                //判断单个bion with ms2 trail匹配是否成功
                //判断连续ion
                //大于某阈值后的连续离子，大于ms2的中平均数？中位数？大于peptide的质量的一半？
                //根据ms2 trail数据的散列度设置阈值？
                if (bpos.size() > bMatchSize) {

                    if (mapCurRTBionMS2Trail == null) {
                        mapCurRTBionMS2Trail = new TreeMap<>();
                    }
                    mapCurRTBionMS2Trail.put((iCharge-1)*len+i +1 ,lMSTwotrails.get(ibposFront));
                    if(Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        bIonMatchTrail[(iCharge-1)*len+i] = lMSTwotrails.get(ibposFront);
                    } else//匹配多个MS2RT
                    {
                        bIonMatchTrail[(iCharge-1)*len+i] = CompateTwoTrail(lMSTwotrails.get(ibposFront), bIonMatchTrail[(iCharge-1)*len+i], b_ionsMZ);
                    }
                }else
                {
                    if(Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        bIonMatchTrail[(iCharge-1)*len+i] = null;
                    }
                }

                if (bIonMatchTrail[(iCharge-1)*len+i] != null)//bpos.size()>bMatchSize)
                {
                    if ((bipos.size() > 0) && (bipos.get(bipos.size() - 1) == (iCharge-1)*len+i - 1)) {
                        iBConsective[iCharge-1]++;
                    }

                    bipos.add((iCharge-1)*len+i);

                    if (mapBionTrail == null) {
                        mapBionTrail = new TreeMap<>();
                    }
                    LinkedList<MSTwoTrail> listMS2Trail = mapBionTrail.get((iCharge-1)*len+i);
                    if (listMS2Trail == null) {
                        listMS2Trail = new LinkedList<MSTwoTrail>();
                    }
                    listMS2Trail.add(bIonMatchTrail[(iCharge-1)*len+i]);
                    mapBionTrail.put((iCharge-1)*len+i, listMS2Trail);

                    dMatch = getdMatchIonScore(bIonMatchTrail[(iCharge-1)*len+i], bIonMaxMass, dMatch, maxPeakArea, iBionCurCount);

                    dMatch = getdMatchMultiIonScore(lMSTwotrails, bipos, bpos, dMatch);

                    //用于输出中间值

                    dlogBionPeakAreaSum[iCharge-1] += Math.log10(bIonMatchTrail[(iCharge-1)*len+i].getPeakArea());
                    dBionMassErrorSum[iCharge-1] += 1.0 - Math.pow((((bIonMatchTrail[(iCharge-1)*len+i].getMzApex() - b_ionsMZ) / (bIonMatchTrail[(iCharge-1)*len+i].getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                    dBionPeakHalfSum[iCharge-1] += Math.max(Math.log10(50.0 / bIonMatchTrail[(iCharge-1)*len+i].getPeakAreaHalfRank()), 0.0);
                    dBionRetentionTimeErrorSum[iCharge-1] += 1 - Math.pow((bIonMatchTrail[(iCharge-1)*len+i].getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
                    if(arrIntensity!=null)
                    {
//                        System.out.println("B arrIntensity "+pep.composition+"  "+msOneTrail.getZ());

                        sumXY += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * arrIntensity[i][iCharge -1 ];
                        sumYY += arrIntensity[i][iCharge -1 ] * arrIntensity[i][iCharge -1 ];
                        sumXX += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea();

                        dbsumXY += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * arrIntensity[i][iCharge -1 ];
                        dbsumYY += arrIntensity[i][iCharge -1 ] * arrIntensity[i][iCharge -1 ];
                        dbsumXX += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea();

                        if(iCharge==1 & i<len*0.6)
                        {
                            db60sumXY += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * arrIntensity[i][iCharge -1 ];
                            db60sumYY += arrIntensity[i][iCharge -1 ] * arrIntensity[i][iCharge -1 ];
                            db60sumXX += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea();

                        }

                        if (arrIntensity[i][iCharge -1 ] > 10 * Double.MIN_VALUE ) {
                            sumXYMatchPredict += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * arrIntensity[i][iCharge -1 ];
                            sumYYMatchPredict += arrIntensity[i][iCharge -1 ] * arrIntensity[i][iCharge -1 ];
                            sumXXMatchPredict += bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea() * bIonMatchTrail[(iCharge - 1) * len + i].getPeakArea();

                        }
                    }
//                    else
//                    {
//                        if(debugOutput++<10) {
//                            System.out.println(" no arrIntensity "+pep.composition+"  "+msOneTrail.getZ());
//                        }
//                    }

                }else
                {

                    if(arrIntensity!=null )
                    {
//                        System.out.println("   bnotin----"+i+"  "+pep.composition+"   "+ msOneTrail.getZ());

//                        System.out.println("arrIntensity:" + JSONWriter.valueToString(arrIntensity));
//                        sumXY += Utils.NoneMatchedIntensity * arrIntensity[i][iCharge -1 ];
//                        sumXX += Utils.NoneMatchedIntensity * Utils.NoneMatchedIntensity;
                        sumYY += arrIntensity[i][iCharge -1] * arrIntensity[i][iCharge -1 ];//做cosine时，没有匹配Y值（不是y离子）也要加

                        dbsumYY += arrIntensity[i][iCharge -1] * arrIntensity[i][iCharge -1 ];//做cosine时，没有匹配Y值（不是y离子）也要加

                        if(iCharge==1 & i<len*0.6)
                        {
                            db60sumYY += arrIntensity[i][iCharge -1] * arrIntensity[i][iCharge -1 ];//做cosine时，没有匹配Y值（不是y离子）也要加

                        }
                    }
                }
                if (ypos.size() > yMatchSize) {

                    if (mapCurRTYionMS2Trail == null) {
                        mapCurRTYionMS2Trail = new TreeMap<>();
                    }
                    mapCurRTYionMS2Trail.put((iCharge-1)*len+i +1 ,lMSTwotrails.get(iyposFront));
                    if(Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        yIonMatchTrail[(iCharge-1)*len+len - i - 1] = lMSTwotrails.get(iyposFront);
                    } else//匹配多个MS2RT
                    {
                        yIonMatchTrail[(iCharge-1)*len+len - i - 1] = CompateTwoTrail(lMSTwotrails.get(iyposFront), yIonMatchTrail[(iCharge-1)*len+len - i - 1], y_ionsMZ);
                    }
                }else
                {
                    if(Utils.bzMatchOneRT)//只匹配一个ms2RT
                    {
                        yIonMatchTrail[(iCharge-1)*len+len - i - 1] = null;
                    }
                }

                if(yIonMatchTrail[(iCharge-1)*len+len - i - 1]!=null)
                {
                    if ((yipos.size() > 0) && (yipos.get(yipos.size() - 1) == (iCharge-1)*len+len - i)) {
                        iYConsective[iCharge-1]++;
                    }


                    yipos.add((iCharge-1)*len+len - i - 1);

                    if (mapYionTrail == null) {
                        mapYionTrail = new TreeMap<>();
                    }
                    LinkedList<MSTwoTrail> listMS2Trail = mapYionTrail.get((iCharge-1)*len+len - i - 1);
                    if (listMS2Trail == null) {
                        listMS2Trail = new LinkedList<MSTwoTrail>();
                    }
                    listMS2Trail.add(yIonMatchTrail[(iCharge-1)*len+len - i - 1]);
                    mapYionTrail.put((iCharge-1)*len+len - i - 1, listMS2Trail);


                    dMatch = getdMatchIonScore(yIonMatchTrail[(iCharge-1)*len+len - i - 1], yIonMaxMass, dMatch, maxPeakArea, iYionCurCount);
                    dMatch = getdMatchMultiIonScore(lMSTwotrails, yipos, ypos, dMatch);

                    dlogYionPeakAreaSum[iCharge-1] += Math.log10(yIonMatchTrail[(iCharge-1)*len+len - i - 1].getPeakArea());
                    dYionMassErrorSum[iCharge-1] += 1.0 - Math.pow((((yIonMatchTrail[(iCharge-1)*len+len - i - 1].getMzApex() - y_ionsMZ) / (yIonMatchTrail[(iCharge-1)*len+len - i - 1].getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                    dYionPeakHalfSum[iCharge-1] += Math.max(Math.log10(50.0 / yIonMatchTrail[(iCharge-1)*len+len - i - 1].getPeakAreaHalfRank()), 0.0);
                    dYionRetentionTimeErrorSum[iCharge-1] += 1 - Math.pow((yIonMatchTrail[(iCharge-1)*len+len - i - 1].getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
                    if(arrIntensity!=null ) {
//                        System.out.println("Y arrIntensity "+pep.composition+"  "+msOneTrail.getZ());

                        sumXY += yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea() * arrIntensity[i][iCharge + 1];
                        sumYY += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                        sumXX += yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea() * yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea();

                        dysumXY += yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea() * arrIntensity[i][iCharge + 1];
                        dysumYY += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                        dysumXX += yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea() * yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea();

                        if (arrIntensity[i][iCharge + 1] > 10 * Double.MIN_VALUE ) {
                            sumXYMatchPredict += yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea() * arrIntensity[i][iCharge + 1];
                            sumYYMatchPredict += arrIntensity[i][iCharge + 1] * arrIntensity[i][iCharge + 1];
                            sumXXMatchPredict += yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea() * yIonMatchTrail[(iCharge - 1) * len + len - i - 1].getPeakArea();


                        }

                    }
//                    else
//                    {
//                        if(debugOutput++<10) {
//                            System.out.println(" no arrIntensity "+pep.composition+"  "+msOneTrail.getZ());
//                        }
//                    }
                }else
                {
                    if(arrIntensity!=null )
                    {
//                        System.out.println("   ynotin----"+i+"  "+pep.composition+"   "+ msOneTrail.getZ());

//                        System.out.println("arrIntensity:" + JSONWriter.valueToString(arrIntensity));
//                        sumXY += Utils.NoneMatchedIntensity * arrIntensity[i][iCharge + 1];
//                        sumXX += Utils.NoneMatchedIntensity * Utils.NoneMatchedIntensity;

                        sumYY += arrIntensity[i][iCharge+ 1] * arrIntensity[i][iCharge+ 1];
                        dysumYY += arrIntensity[i][iCharge+ 1] * arrIntensity[i][iCharge+ 1];
                    }

                }

            }

//            if(msOneTrail.getId()==24474)
//                System.out.println("ok");
            bSum[iCharge-1] = bipos.size();
            ySum[iCharge-1] = yipos.size();
            bySum = bySum + bSum[iCharge-1]+ySum[iCharge-1];

            if (FastMath.abs(dbsumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dbsumYY) < 10 * Double.MIN_VALUE) {
                dBionPearsonCorrelationWithPredictMSMS[iCharge-1] =  0.0;
            }else
            {
                dBionPearsonCorrelationWithPredictMSMS[iCharge-1] = dbsumXY/FastMath.sqrt(dbsumXX*dbsumYY);
            }

            if (FastMath.abs(dysumXX) < 10 * Double.MIN_VALUE || FastMath.abs(dysumYY) < 10 * Double.MIN_VALUE) {
                dYionPearsonCorrelationWithPredictMSMS[iCharge-1] =  0.0;
            }else
            {
                dYionPearsonCorrelationWithPredictMSMS[iCharge-1] = dysumXY/FastMath.sqrt(dysumXX*dysumYY);
            }

            if (bySum > Utils.thresholdFragmentCount/*0*/) {
                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance) {
                    dMatch += dMatch * log10(bySum) / log10(lMSTwotrails.size());
//                dMatch = dMatch * bySum / lMSTwotrails.size();

                }
//                strb += "|B"+iCharge+":" + bipos.stream().map(x -> x + 1).collect(Collectors.toList()).toString();
//                stry += "|Y"+iCharge+":" + yipos.stream().map(x -> len - x).collect(Collectors.toList()).toString();

//            dMatch = dMatch * log10(bySum) / log10(2*pep.b_ions.length);
//            dMatch = dMatch * log10(bySum);


            } else {
                dMatch = 0.0;//未匹配到
            }
        }

        if (bySum >  Utils.thresholdFragmentCount/*0*/) {
            String strDecoy = "target";
//                if(pep.id.contains("Decoy"))//DeBruijnsp
            if (pep.id.contains("DeBruijn"))//>DeBruijnsp
            {
                strDecoy = "decoy";
            }
            String mstwortapex = "" + lMSTwotrails.get(0).getRtApex();
            matchPep = new MatchPeptideWithMultiCharge();
            matchPep.setValues(
                    dMatch, pep.id, strDecoy, pep.composition, pep.mass, msOneTrail,
                    mstwortapex,  bSum, ySum,  dMatch, strb
                    , stry, bipos, yipos);
            matchPep.setPeakAreaMSErrorPeakHalfRankWithMultiCharge(dlogBionPeakAreaSum, dlogYionPeakAreaSum,
                    dBionMassErrorSum, dYionMassErrorSum, dBionPeakHalfSum, dYionPeakHalfSum,
                    dBionRetentionTimeErrorSum, dYionRetentionTimeErrorSum, pep.dMutationRate);

            matchPep.setConsetivePeakWithMultiCharge(iBConsective, iYConsective);
            matchPep.setPep(pep);

            matchPep.setMapBionMS2Trail(mapBionTrail);
            matchPep.setMapYionMS2Trail(mapYionTrail);
            matchPep.mapCurRTBionMS2Trail = mapCurRTBionMS2Trail;
            matchPep.mapCurRTYionMS2Trail = mapCurRTYionMS2Trail;

            matchPep.arr_dBionPearsonCorrelationWithPredictMSMS = dBionPearsonCorrelationWithPredictMSMS;
            matchPep.arr_dYionPearsonCorrelationWithPredictMSMS = dYionPearsonCorrelationWithPredictMSMS;


            //B 1 ION 60%
            if (FastMath.abs(db60sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(db60sumYY) < 10 * Double.MIN_VALUE) {
                matchPep.dBion60PearsonCorrelationWithPredictMSMS= 0.0;
            }else
            {
                matchPep.dBion60PearsonCorrelationWithPredictMSMS= db60sumXY/FastMath.sqrt(db60sumXX*db60sumYY);


            }

            //matched and predict cos
            if (FastMath.abs(sumXXMatchPredict) < 10 * Double.MIN_VALUE || FastMath.abs(sumYYMatchPredict) < 10 * Double.MIN_VALUE) {
                matchPep.dCosWithPredictMSMSMatchedPredict= 0.0;
            }else
            {
                matchPep.dCosWithPredictMSMSMatchedPredict=  sumXYMatchPredict/FastMath.sqrt(sumXXMatchPredict*sumYYMatchPredict);


            }

            if (FastMath.abs(sumXX) < 10 * Double.MIN_VALUE || FastMath.abs(sumYY) < 10 * Double.MIN_VALUE) {
                matchPep.dPearsonCorrelationWithPredictMSMS= 0.0;
            }else
            {
                matchPep.dPearsonCorrelationWithPredictMSMS= sumXY/FastMath.sqrt(sumXX*sumYY);
//                if(debugOutput++<100) {
//                    System.out.println("In matchPepCos:" + matchPep.dPearsonCorrelationWithPredictMSMS);
//                }

            }


        }
        //       return dMatch;
        return  matchPep ;

    }

/*
    public  static  MatchPeptdide PepMatchIonCountCrossRT(Peptide pep, List<MSTwoTrail> lMSTwotrails, double maxPeakArea,
                                                   Enums.ScoreMode scoreMode, BufferedWriter writer, int icount,
                                                   double ms1rt, double ms1Mass, long ms1id,double ms1qualityScore,
                                                    MSTwoTrail[] bIonMatchTrail,MSTwoTrail[] yIonMatchTrail) {

        double dMatch = 1.0;
        double dlogBionPeakAreaSum = 0.0;
        double dlogYionPeakAreaSum = 0.0;

        double dBionMassErrorSum = 0.0;
        double dYionMassErrorSum = 0.0;
        double dBionRetentionTimeErrorSum = 0.0;
        double dYionRetentionTimeErrorSum = 0.0;

        double dBionPeakHalfSum = 0.0;
        double dYionPeakHalfSum = 0.0;

        int iBConsective = 0;
        int iYConsective = 0;

        TreeMap<Integer, LinkedList<MSTwoTrail>>  mapBionTrail = null;
        TreeMap<Integer, LinkedList<MSTwoTrail>>  mapYionTrail = null;



        if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances || Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
        {
            dMatch = 0.0;//add all abundance
        }
        int len = pep.b_ions.length;
//        int iCurrentPos = 0;
        bpos.clear();
        ypos.clear();

        bipos.clear();
        yipos.clear();


        int ibposFront = 0;
        int ibposBack = 0;
        int iyposFront = 0;
        int iyposBack = 0;

        int ibnearpos = 0;
        int iynearpos = 0;

        int bMatchSize;
        int yMatchSize;

//        if(ms1id ==6671)
//            System.out.println("test");
        for(int i = 0; i < len; ++i) {

            bMatchSize = bpos.size();
            yMatchSize = ypos.size();

            double minBIonDistance = Utils.thresholdPPM;
            double minYIonDistance = Utils.thresholdPPM;

            double bIonMinusMass = pep.b_ions[i]*(1 - Utils.thresholdPPM);
            double bIonMaxMass = pep.b_ions[i]*(1 + Utils.thresholdPPM);

            int  iBionCurCount = 0;
            int  iYionCurCount = 0;



            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( bIonMaxMass > (lMSTwotrails.get(0).getMzApex()) )
                    && ( bIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
//                if(i>0)
//                {
//                    if(bIonMinusMass<pep.y_ions[len-i])//根据上一次的找的y ion 位置，更新b ion 需要找的空间
//                    {
//                        ibposFront = binarySearch0(lMSTwotrails,ibposFront,iyposFront,bIonMinusMass);
//                    }
//                    else
//                    {
//                        ibposFront = binarySearch0(lMSTwotrails,iyposFront,lMSTwotrails.size(),bIonMinusMass);
//                    }
//
//
//                }else
//                {
                ibposFront = binarySearch0(lMSTwotrails,ibposFront,lMSTwotrails.size(),bIonMinusMass);//全新查找
//                }

//            ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size()-1,pep.b_ions[i]);
                ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposFront);

                ///更新bpos里面的内容
                iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass,bpos);

//                ibposBack = binarySearch0(lMSTwotrails,ibnearpos,lMSTwotrails.size(),bIonMaxMass);
//                if (ibposBack < 0)//在区间内
//                {
//                    ibposBack = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposBack);
//                    if (ibposBack>0)
//                        ibposBack--;
//                }
//
//                iBionCurCount = ibposBack - ibnearpos + 1;
//                if( iBionCurCount>0 )
//                {
//                    bpos.add(ibnearpos);
//                }
                if (bpos.size() > 0)
                    ibposFront = bpos.get(bpos.size() - 1);
                if (ibposFront < 0)
                    ibposFront = ibnearpos;
            }


            double yIonMinusMass = pep.y_ions[len-i-1]*(1 - Utils.thresholdPPM);
            double yIonMaxMass = pep.y_ions[len-i-1]*(1 + Utils.thresholdPPM);


            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( yIonMaxMass > (lMSTwotrails.get(0).getMzApex()) )
                    && ( yIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
                //根据前面比较的结果，search 对应的y离子
//                if (bIonMinusMass > yIonMinusMass)//根据以前的找的b ion 位置，更新y ion 需要找的空间
//                {
//                    iyposFront = binarySearch0(lMSTwotrails,iyposFront,ibposFront,yIonMinusMass);
//
//                }else
//                {
                iyposFront = binarySearch0(lMSTwotrails,iyposFront,lMSTwotrails.size(),yIonMinusMass);
//                }
                iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);

                iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass,ypos);

//                int iyposnext = iynearpos;
//                //若是在精度范围区间，开始计数
//                for(  ;iyposnext<lMSTwotrails.size() &&
//                        lMSTwotrails.get(iyposnext).getMzApex() <= yIonMaxMass;
//                      iyposnext++)
//                {
//                    iYionCurCount++;
//                }
//                if (iYionCurCount > 0) {
//                    ypos.add(iynearpos);
//                }

//                iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);
//                iyposBack = binarySearch0(lMSTwotrails, iynearpos, lMSTwotrails.size(), yIonMaxMass);
//                if (iyposBack < 0)//在区间内
//                {
//                    iyposBack = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposBack);
//                    if (iyposBack > 0)
//                        iyposBack--;
//                }

//                iYionCurCount = iyposBack - iynearpos  + 1;
//                if (iYionCurCount > 0) {
//                    ypos.add(iynearpos);
//                }


                if(ypos.size()>0)
                    iyposFront = ypos.get(ypos.size()-1);
                if(iyposFront < 0)
                    iyposFront = iynearpos;
            }


            //判断单个bion with ms2 trail匹配是否成功
            //判断连续ion
            //大于某阈值后的连续离子，大于ms2的中平均数？中位数？大于peptide的质量的一半？
            //根据ms2 trail数据的散列度设置阈值？
            //
//            decoy	LTCNYLVTPFEGAFYNFSLNLK	|B:[2, 4, 8, 12, 13, 14, 15, 20]	|Y:[20, 19, 18, 17, 9, 8, 7]
//            decoy	VICQEPHISGALTSAAEEAR	|B:[3, 4, 5, 9, 14, 16]	|Y:[18, 17, 16, 15, 11, 9, 3]
//            decoy	FAVIDSFTRPDQVQHTWTK	|B:[2, 3, 6, 9, 10, 12, 13, 15]	|Y:[17, 16, 11, 8, 4, 3]
//            decoy	LATSSLDHAGGATGAAGAAGAAGGGGAAELR	|B:[7, 11, 13, 15, 16, 21, 29]	|Y:[29, 28, 16, 13, 12, 10, 9, 7, 3]
//            decoy	NIGLILK	|B:[2, 3, 4, 5, 6]	|Y:[5, 4, 3, 2, 1]
//            decoy	NIGLILK	|B:[2, 3, 4, 5, 6]	|Y:[5, 4, 3, 2, 1]
//            decoy	AEIFHLEQAVHCATLPNDGTIMVR	|B:[2, 3, 4, 6, 7, 11, 12, 13, 14, 17, 18]	|Y:[16, 11, 8]
//            decoy	AEIFHLEQAVHCATLPNDGTIMVR	|B:[2, 3, 4, 6, 7, 11, 12, 13, 14, 17, 18]	|Y:[16, 11, 8]
//            decoy	FPSAEFICPISVTSGFAAFK	|B:[2, 7, 10, 12, 13, 18]	|Y:[18, 17, 15, 7, 6, 3, 2, 1]
//            decoy	FPSAEFICPISVTSGFAAFK	|B:[2, 7, 10, 12, 13, 18]	|Y:[18, 17, 15, 7, 6, 3, 2, 1]
//            decoy	HLAAEVER	|B:[2, 3, 4, 7]	|Y:[6, 5, 4, 3, 2, 1]
//            decoy	HLDIEVTQQGAPINPDGTK	|B:[2, 4, 5, 6, 7, 10, 11, 14]	|Y:[16, 11, 10, 8, 7]
//            decoy	SLNSFYLGSVASAPAMALSETR	|B:[5, 8, 13, 14, 16, 17, 18]	|Y:[17, 15, 14, 13, 11, 8, 7]
//            decoy	SGLPSEAEISVSLASYMEVDK	|B:[3, 5, 7, 9, 11, 15, 16, 17]	|Y:[19, 10, 7, 5, 4, 2]
//            decoy	TIGATQVYTANQSVAGGAQAPMK	|B:[3, 9, 12, 14, 15, 16, 17, 19]	|Y:[20, 15, 12, 11, 9, 8]
//            decoy	VALFSESPDPFVIEDYSEVVVK	|B:[3, 5, 6, 7, 10, 12, 13, 15]	|Y:[20, 19, 14, 12, 11, 10]
//            decoy	ISIPYLTNHPDGNGMNITGK	|B:[2, 3, 5, 12, 14, 17, 18]	|Y:[18, 17, 14, 13, 11, 5, 2]
//            decoy	NDLGLLEAR	|B:[2, 3, 4, 5]	|Y:[6, 5, 3, 2, 1]
//            decoy	ITDFHGDVLELIDGETSHGFLK	|B:[2, 3, 9, 10, 11, 13, 15]	|Y:[20, 19, 18, 16, 14, 13, 9]
//            decoy	VQFTDASVLLAQELR	|B:[2, 3, 6, 8, 9, 11]	|Y:[13, 11, 10, 9, 7, 3]
//            decoy	LDIEVQGAQNEIMAHQQK	|B:[2, 6, 9, 12, 13, 14, 15]	|Y:[16, 15, 10, 9, 8, 5]
             //
            if( bpos.size()>bMatchSize)
            {
                bIonMatchTrail[i] = CompateTwoTrail(lMSTwotrails.get(ibposFront),bIonMatchTrail[i], pep.b_ions[i]);
            }

            if( bIonMatchTrail[i]!=null)//bpos.size()>bMatchSize)
            {
                if ((bipos.size()>0) && (bipos.get(bipos.size()-1)==i-1))
                {
                    iBConsective++;
                }

                bipos.add(i);

                if(mapBionTrail==null)
                {
                    mapBionTrail = new TreeMap<>();
                }
                LinkedList<MSTwoTrail> listMS2Trail =  mapBionTrail.get(i);
                if(listMS2Trail==null)
                {
                    listMS2Trail = new LinkedList<MSTwoTrail>();
                }
                listMS2Trail.add(bIonMatchTrail[i]);
                mapBionTrail.put(i,listMS2Trail);

                dMatch = getdMatchIonScore(bIonMatchTrail[i], bIonMaxMass, dMatch, maxPeakArea, iBionCurCount);

                dMatch = getdMatchMultiIonScore(lMSTwotrails, bipos,bpos, dMatch);

                //用于输出中间值

                dlogBionPeakAreaSum += Math.log10(bIonMatchTrail[i].getPeakArea());
                dBionMassErrorSum += 1.0-Math.pow((((bIonMatchTrail[i].getMzApex() - pep.b_ions[i])/(bIonMatchTrail[i].getMzApex() / 1000000.0))/(Utils.thresholdPPM*1000000.0)),2);
                dBionPeakHalfSum += Math.max(Math.log10(50.0/bIonMatchTrail[i].getPeakAreaHalfRank()),0.0);
                dBionRetentionTimeErrorSum += 1-Math.pow((bIonMatchTrail[i].getRtApex() - ms1rt)/(Utils.thresholdRT),2);

            }
            if( ypos.size()>yMatchSize)
            {
                yIonMatchTrail[len-i-1] = CompateTwoTrail(lMSTwotrails.get(iyposFront),yIonMatchTrail[len-i-1], pep.y_ions[len-i-1]);
            }
            if(yIonMatchTrail[len-i-1]!=null)
            {
                if ((yipos.size()>0) && (yipos.get(yipos.size()-1)==len-i))
                {
                    iYConsective++;
                }


                yipos.add(len-i-1);

                if(mapYionTrail==null)
                {
                    mapYionTrail = new TreeMap<>();
                }
                LinkedList<MSTwoTrail> listMS2Trail =  mapYionTrail.get(len-i-1);
                if(listMS2Trail==null)
                {
                    listMS2Trail = new LinkedList<MSTwoTrail>();
                }
                listMS2Trail.add(yIonMatchTrail[len-i-1] );
                mapYionTrail.put(len-i-1,listMS2Trail);


                dMatch = getdMatchIonScore( yIonMatchTrail[len-i-1], yIonMaxMass, dMatch, maxPeakArea, iYionCurCount);
                dMatch = getdMatchMultiIonScore(lMSTwotrails, yipos,ypos, dMatch);

                dlogYionPeakAreaSum += Math.log10(yIonMatchTrail[len-i-1].getPeakArea());
                dYionMassErrorSum += 1.0-Math.pow((((yIonMatchTrail[len-i-1].getMzApex() -  pep.y_ions[len-i-1])/(yIonMatchTrail[len-i-1].getMzApex() / 1000000.0))/(Utils.thresholdPPM*1000000.0)),2);
                dYionPeakHalfSum  += Math.max(Math.log10(50.0/yIonMatchTrail[len-i-1].getPeakAreaHalfRank()),0.0);
                dYionRetentionTimeErrorSum += 1-Math.pow((yIonMatchTrail[len-i-1].getRtApex() - ms1rt)/(Utils.thresholdRT),2);

//                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
//                {
//                    //                    dMatch += Math.log10(lMSTwotrails.get(iyposFront).getPeakArea())/iYionCurCount;//(Math.log10(maxPeakArea)*iYionCurCount);//得到对应y ion位置的峰面积
//                  dMatch += Math.log10(lMSTwotrails.get(iyposFront).getPeakArea());//得到对应b ion位置的峰
//
//                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
//                {
//                    dMatch*=(1-minYIonDistance);//得到理论mz与实际质荷比距离
//
//                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
//                {
//                    dMatch += (1-minYIonDistance)*Math.log10(lMSTwotrails.get(iyposFront).getPeakArea())/iYionCurCount;//(Math.log10(maxPeakArea)*iYionCurCount);//得到理论mz与实际质荷比距离，对应位置的峰面积
//
//                }
            }
        }
        MatchPeptdide matchPep = null;
        int bySum = bpos.size() + ypos.size();
        if(bySum>0)
        {
            if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance) {
                dMatch = dMatch * log10(bySum) / log10(lMSTwotrails.size());
//                dMatch = dMatch * bySum / lMSTwotrails.size();

            }
//            dMatch = dMatch * log10(bySum) / log10(2*pep.b_ions.length);
//            dMatch = dMatch * log10(bySum);
            if(dMatch>0.0) {
                String strDecoy="target";
//                if(pep.id.contains("Decoy"))//DeBruijnsp
                if(pep.id.contains("DeBruijnsp"))//>DeBruijnsp
                {
                    strDecoy = "decoy";
                }
                String mstwortapex =""+lMSTwotrails.get(0).getRtApex();
                String strb= "|B:" + bipos.stream().map(x->x+1).collect(Collectors.toList()).toString();
                String stry= "|Y:" + yipos.stream().map(x->len-x).collect(Collectors.toList()).toString();
                matchPep = new MatchPeptdide();
                matchPep.setValues(
                        ms1id, dMatch,pep.id,strDecoy,pep.composition,pep.mass,ms1Mass,ms1rt,ms1qualityScore,
                        mstwortapex,"" + bySum,""+dMatch,strb
                        ,stry,bipos,yipos);
                matchPep.setPeakAreaMSErrorPeakHalfRank(dlogBionPeakAreaSum,dlogYionPeakAreaSum,
                        dBionMassErrorSum,dYionMassErrorSum,dBionPeakHalfSum,dYionPeakHalfSum,
                        dBionRetentionTimeErrorSum,dYionRetentionTimeErrorSum,pep.dMutationRate);

                matchPep.setConsetivePeak(iBConsective,iYConsective);
                matchPep.setPep(pep);





//                for(Integer ib:bpos)
//                {
//                    lMSTwotrails.get(ib).addMatchPeptide(matchPep);
//                }
//                for(Integer iy:ypos)
//                {
//                    lMSTwotrails.get(iy).addMatchPeptide(matchPep);
//                }

//                int iL = 0;
//                for(Integer ibp:bipos)
//                {
//                    matchPep.addmapBionMS2Trail(ibp,lMSTwotrails.get(bpos.get(iL++)));
//                }
//                iL = 0;
//                for(Integer iyp:yipos)
//                {
//                    matchPep.addmapYionMS2Trail(iyp,lMSTwotrails.get(ypos.get(iL++)));
//                }
                matchPep.setMapBionMS2Trail(mapBionTrail);
                matchPep.setMapYionMS2Trail(mapYionTrail);


            }

        }
        else
            dMatch = 0.0;//未匹配到

        //       return dMatch;
        return  matchPep ;

    }
*/
    /*
    public  static  MatchPeptdide PepMatchIonCount(Peptide pep, List<MSTwoTrail> lMSTwotrails, double maxPeakArea,
                                                   Enums.ScoreMode scoreMode, BufferedWriter writer, int icount,
                                                   double ms1rt, double ms1Mass, long ms1id,double ms1qualityScore) {

//        public DbBaseMatch PepMatch(Peptide pep, List<MSTwoTrail> lMSTwotrails,double maxPeakArea) {
        double dMatch = 1.0;
        double dlogBionPeakAreaSum = 0.0;
        double dlogYionPeakAreaSum = 0.0;

        double dBionMassErrorSum = 0.0;
        double dYionMassErrorSum = 0.0;

        double dBionPeakHalfSum = 0.0;
        double dYionPeakHalfSum = 0.0;

        double dBionRetentionTimeErrorSum = 0.0;
        double dYionRetentionTimeErrorSum = 0.0;

        int iBConsective = 0;
        int iYConsective = 0;

        if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances || Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
        {
            dMatch = 0.0;//add all abundance
        }
        int len = pep.b_ions.length;
//        int iCurrentPos = 0;
        bpos.clear();
        ypos.clear();

        bipos.clear();
        yipos.clear();


        int ibposFront = 0;
        int ibposBack = 0;
        int iyposFront = 0;
        int iyposBack = 0;

        int ibnearpos = 0;
        int iynearpos = 0;

        int bMatchSize;
        int yMatchSize;

//        if(ms1id ==6671)
//            System.out.println("test");
        for(int i = 0; i < len; ++i) {

            bMatchSize = bpos.size();
            yMatchSize = ypos.size();

            double minBIonDistance = Utils.thresholdPPM;
            double minYIonDistance = Utils.thresholdPPM;

            double bIonMinusMass = pep.b_ions[i]*(1 - Utils.thresholdPPM);
            double bIonMaxMass = pep.b_ions[i]*(1 + Utils.thresholdPPM);

            int  iBionCurCount = 0;
            int  iYionCurCount = 0;



            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( bIonMaxMass > (lMSTwotrails.get(0).getMzApex()) )
                    && ( bIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
//                if(i>0)
//                {
//                    if(bIonMinusMass<pep.y_ions[len-i])//根据上一次的找的y ion 位置，更新b ion 需要找的空间
//                    {
//                        ibposFront = binarySearch0(lMSTwotrails,ibposFront,iyposFront,bIonMinusMass);
//                    }
//                    else
//                    {
//                        ibposFront = binarySearch0(lMSTwotrails,iyposFront,lMSTwotrails.size(),bIonMinusMass);
//                    }
//
//
//                }else
//                {
                    ibposFront = binarySearch0(lMSTwotrails,ibposFront,lMSTwotrails.size(),bIonMinusMass);//全新查找
//                }

//            ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size()-1,pep.b_ions[i]);
                ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposFront);

                ///更新bpos里面的内容
                iBionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, ibnearpos, bIonMaxMass,bpos);

//                ibposBack = binarySearch0(lMSTwotrails,ibnearpos,lMSTwotrails.size(),bIonMaxMass);
//                if (ibposBack < 0)//在区间内
//                {
//                    ibposBack = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposBack);
//                    if (ibposBack>0)
//                        ibposBack--;
//                }
//
//                iBionCurCount = ibposBack - ibnearpos + 1;
//                if( iBionCurCount>0 )
//                {
//                    bpos.add(ibnearpos);
//                }
                if (bpos.size() > 0)
                    ibposFront = bpos.get(bpos.size() - 1);
                if (ibposFront < 0)
                    ibposFront = ibnearpos;
            }


            double yIonMinusMass = pep.y_ions[len-i-1]*(1 - Utils.thresholdPPM);
            double yIonMaxMass = pep.y_ions[len-i-1]*(1 + Utils.thresholdPPM);


            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( yIonMaxMass > (lMSTwotrails.get(0).getMzApex()) )
                    && ( yIonMinusMass < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
                //根据前面比较的结果，search 对应的y离子
//                if (bIonMinusMass > yIonMinusMass)//根据以前的找的b ion 位置，更新y ion 需要找的空间
//                {
//                    iyposFront = binarySearch0(lMSTwotrails,iyposFront,ibposFront,yIonMinusMass);
//
//                }else
//                {
                iyposFront = binarySearch0(lMSTwotrails,iyposFront,lMSTwotrails.size(),yIonMinusMass);
//                }
                iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);

                iYionCurCount = getiIonCurCountWithCheckIsotope(lMSTwotrails, iynearpos, yIonMaxMass,ypos);

//                int iyposnext = iynearpos;
//                //若是在精度范围区间，开始计数
//                for(  ;iyposnext<lMSTwotrails.size() &&
//                        lMSTwotrails.get(iyposnext).getMzApex() <= yIonMaxMass;
//                      iyposnext++)
//                {
//                    iYionCurCount++;
//                }
//                if (iYionCurCount > 0) {
//                    ypos.add(iynearpos);
//                }

//                iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);
//                iyposBack = binarySearch0(lMSTwotrails, iynearpos, lMSTwotrails.size(), yIonMaxMass);
//                if (iyposBack < 0)//在区间内
//                {
//                    iyposBack = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposBack);
//                    if (iyposBack > 0)
//                        iyposBack--;
//                }

//                iYionCurCount = iyposBack - iynearpos  + 1;
//                if (iYionCurCount > 0) {
//                    ypos.add(iynearpos);
//                }


                if(ypos.size()>0)
                    iyposFront = ypos.get(ypos.size()-1);
                if(iyposFront < 0)
                    iyposFront = iynearpos;
            }


            //判断单个bion with ms2 trail匹配是否成功
            //判断连续ion
            //大于某阈值后的连续离子，大于ms2的中平均数？中位数？大于peptide的质量的一半？
            //根据ms2 trail数据的散列度设置阈值？
            //
//            decoy	LTCNYLVTPFEGAFYNFSLNLK	|B:[2, 4, 8, 12, 13, 14, 15, 20]	|Y:[20, 19, 18, 17, 9, 8, 7]
//            decoy	VICQEPHISGALTSAAEEAR	|B:[3, 4, 5, 9, 14, 16]	|Y:[18, 17, 16, 15, 11, 9, 3]
//            decoy	FAVIDSFTRPDQVQHTWTK	|B:[2, 3, 6, 9, 10, 12, 13, 15]	|Y:[17, 16, 11, 8, 4, 3]
//            decoy	LATSSLDHAGGATGAAGAAGAAGGGGAAELR	|B:[7, 11, 13, 15, 16, 21, 29]	|Y:[29, 28, 16, 13, 12, 10, 9, 7, 3]
//            decoy	NIGLILK	|B:[2, 3, 4, 5, 6]	|Y:[5, 4, 3, 2, 1]
//            decoy	NIGLILK	|B:[2, 3, 4, 5, 6]	|Y:[5, 4, 3, 2, 1]
//            decoy	AEIFHLEQAVHCATLPNDGTIMVR	|B:[2, 3, 4, 6, 7, 11, 12, 13, 14, 17, 18]	|Y:[16, 11, 8]
//            decoy	AEIFHLEQAVHCATLPNDGTIMVR	|B:[2, 3, 4, 6, 7, 11, 12, 13, 14, 17, 18]	|Y:[16, 11, 8]
//            decoy	FPSAEFICPISVTSGFAAFK	|B:[2, 7, 10, 12, 13, 18]	|Y:[18, 17, 15, 7, 6, 3, 2, 1]
//            decoy	FPSAEFICPISVTSGFAAFK	|B:[2, 7, 10, 12, 13, 18]	|Y:[18, 17, 15, 7, 6, 3, 2, 1]
//            decoy	HLAAEVER	|B:[2, 3, 4, 7]	|Y:[6, 5, 4, 3, 2, 1]
//            decoy	HLDIEVTQQGAPINPDGTK	|B:[2, 4, 5, 6, 7, 10, 11, 14]	|Y:[16, 11, 10, 8, 7]
//            decoy	SLNSFYLGSVASAPAMALSETR	|B:[5, 8, 13, 14, 16, 17, 18]	|Y:[17, 15, 14, 13, 11, 8, 7]
//            decoy	SGLPSEAEISVSLASYMEVDK	|B:[3, 5, 7, 9, 11, 15, 16, 17]	|Y:[19, 10, 7, 5, 4, 2]
//            decoy	TIGATQVYTANQSVAGGAQAPMK	|B:[3, 9, 12, 14, 15, 16, 17, 19]	|Y:[20, 15, 12, 11, 9, 8]
//            decoy	VALFSESPDPFVIEDYSEVVVK	|B:[3, 5, 6, 7, 10, 12, 13, 15]	|Y:[20, 19, 14, 12, 11, 10]
//            decoy	ISIPYLTNHPDGNGMNITGK	|B:[2, 3, 5, 12, 14, 17, 18]	|Y:[18, 17, 14, 13, 11, 5, 2]
//            decoy	NDLGLLEAR	|B:[2, 3, 4, 5]	|Y:[6, 5, 3, 2, 1]
//            decoy	ITDFHGDVLELIDGETSHGFLK	|B:[2, 3, 9, 10, 11, 13, 15]	|Y:[20, 19, 18, 16, 14, 13, 9]
//            decoy	VQFTDASVLLAQELR	|B:[2, 3, 6, 8, 9, 11]	|Y:[13, 11, 10, 9, 7, 3]
//            decoy	LDIEVQGAQNEIMAHQQK	|B:[2, 6, 9, 12, 13, 14, 15]	|Y:[16, 15, 10, 9, 8, 5]
             //

            if(bpos.size()>bMatchSize)
            {
                if ((bipos.size()>0) && (bipos.get(bipos.size()-1)==i-1))
                {
                    iBConsective++;
                }

                bipos.add(i);
                dMatch = getdMatchIonScore(lMSTwotrails, bIonMaxMass, dMatch, ibposFront, maxPeakArea, iBionCurCount);

                dMatch = getdMatchMultiIonScore(lMSTwotrails, bipos,bpos, dMatch);

                //用于输出中间值

                dlogBionPeakAreaSum += Math.log10(lMSTwotrails.get(ibposFront).getPeakArea());
                dBionMassErrorSum += 1.0-Math.pow((((lMSTwotrails.get(ibposFront).getMzApex() -  pep.b_ions[i])/(lMSTwotrails.get(ibposFront).getMzApex() / 1000000.0))/(Utils.thresholdPPM*1000000.0)),2);
                dBionPeakHalfSum += Math.max(Math.log10(50.0/lMSTwotrails.get(ibposFront).getPeakAreaHalfRank()),0.0);
                dBionRetentionTimeErrorSum += 1-Math.pow((lMSTwotrails.get(ibposFront).getRtApex() - ms1rt)/(Utils.thresholdRT),2);


            }
            if(ypos.size()>yMatchSize)
            {
                if ((yipos.size()>0) && (yipos.get(yipos.size()-1)==len-i))
                {
                    iYConsective++;
                }

                yipos.add(len-i-1);
                dMatch = getdMatchIonScore(lMSTwotrails, yIonMaxMass, dMatch, iyposFront, maxPeakArea, iYionCurCount);
                dMatch = getdMatchMultiIonScore(lMSTwotrails, yipos,ypos, dMatch);

                dlogYionPeakAreaSum += Math.log10(lMSTwotrails.get(iyposFront).getPeakArea());
                dYionMassErrorSum += 1.0-Math.pow((((lMSTwotrails.get(iyposFront).getMzApex() -  pep.y_ions[len-i-1])/(lMSTwotrails.get(iyposFront).getMzApex() / 1000000.0))/(Utils.thresholdPPM*1000000.0)),2);
                dYionPeakHalfSum  += Math.max(Math.log10(50.0/lMSTwotrails.get(iyposFront).getPeakAreaHalfRank()),0.0);
                dYionRetentionTimeErrorSum += 1-Math.pow((lMSTwotrails.get(iyposFront).getRtApex() - ms1rt)/(Utils.thresholdRT),2);

//                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
//                {
//                    //                    dMatch += Math.log10(lMSTwotrails.get(iyposFront).getPeakArea())/iYionCurCount;//(Math.log10(maxPeakArea)*iYionCurCount);//得到对应y ion位置的峰面积
//                  dMatch += Math.log10(lMSTwotrails.get(iyposFront).getPeakArea());//得到对应b ion位置的峰
//
//                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
//                {
//                    dMatch*=(1-minYIonDistance);//得到理论mz与实际质荷比距离
//
//                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
//                {
//                    dMatch += (1-minYIonDistance)*Math.log10(lMSTwotrails.get(iyposFront).getPeakArea())/iYionCurCount;//(Math.log10(maxPeakArea)*iYionCurCount);//得到理论mz与实际质荷比距离，对应位置的峰面积
//
//                }
            }
        }
        MatchPeptdide matchPep = null;
        int bySum = bpos.size() + ypos.size();
        if(bySum>0)
        {
            if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance) {
                dMatch = dMatch * log10(bySum) / log10(lMSTwotrails.size());
//                dMatch = dMatch * bySum / lMSTwotrails.size();

            }
//            dMatch = dMatch * log10(bySum) / log10(2*pep.b_ions.length);
//            dMatch = dMatch * log10(bySum);
            if(dMatch>0.0) {
                String strDecoy="target";
//                if(pep.id.contains("Decoy"))//DeBruijnsp
                if(pep.id.contains("DeBruijnsp"))//>DeBruijnsp
                {
                    strDecoy = "decoy";
                }
                String mstwortapex =""+lMSTwotrails.get(0).getRtApex();
                String strb= "|B:" + bipos.stream().map(x->x+1).collect(Collectors.toList()).toString();
                String stry= "|Y:" + yipos.stream().map(x->x+1).collect(Collectors.toList()).toString();
                matchPep = new MatchPeptdide();
                matchPep.setValues(
                        ms1id, dMatch,pep.id,strDecoy,pep.composition,pep.mass,ms1Mass,ms1rt,ms1qualityScore,
                        mstwortapex,"" + bySum,""+dMatch,strb
                        ,stry,bipos,yipos);
                matchPep.setPeakAreaMSErrorPeakHalfRank(dlogBionPeakAreaSum,dlogYionPeakAreaSum,
                        dBionMassErrorSum,dYionMassErrorSum,dBionPeakHalfSum,dYionPeakHalfSum,
                        dBionRetentionTimeErrorSum,dYionRetentionTimeErrorSum,pep.dMutationRate);
                matchPep.setConsetivePeak(iBConsective,iYConsective);

            }

        }
        else
            dMatch = 0.0;//未匹配到

        //       return dMatch;
        return  matchPep ;

    }

     */

    private static double getdMatchMultiIonScore(List<MSTwoTrail> lMSTwotrails, List<Integer> ipos,List<Integer> pos,double dMatch) {
        if(Utils.MutiIonMatchMode== Utils.MutiIonMatchEnum.Consective)
        {
            if(ipos.size() >2 )
            {
                if(Math.abs(ipos.get(ipos.size()-1) - ipos.get(ipos.size()-2))<2)//与前面一个匹配的离子相连吗？
                {
                    dMatch *=1.1;
                }
            }
        }else if(Utils.MutiIonMatchMode== Utils.MutiIonMatchEnum.BackLargeMassSearch)//假设后面的离子匹配概率低一些
        {
  //          dMatch *=Math.sqrt((pos.get(pos.size()-1)+1.0)/ lMSTwotrails.size());//become a double//这里以后优化，质量大一些匹配概率低一些，不是pos里面的位置

        }else if(Utils.MutiIonMatchMode== Utils.MutiIonMatchEnum.ConbineConsectiveAndBackLargeMassSearch)
        {
            if(ipos.size() >2)
            {
                if(Math.abs(ipos.get(ipos.size()-1) - ipos.get(ipos.size()-2))<2)//与前面一个匹配的离子相连吗？
                {
                    dMatch *=1.1;
                }
            }
   //         dMatch *=Math.sqrt((pos.get(pos.size()-1)+1.0)/ lMSTwotrails.size());//become a double

        }
        return dMatch;
    }

    private static double getdMatchIonScore(MSTwoTrail lMSTwotrails, double maxIonMaxMass, double dMatch,  double maxPeakArea, int iIonCurCount) {
        if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
        {
            //                   dMatch += Math.log10(lMSTwotrails.get(ibposFront).getPeakArea())/(Math.log10(maxPeakArea)*iBionCurCount);//得到对应b ion位置的峰面积 iBionCurCount;//
//            dMatch += Math.log10(lMSTwotrails.get(iposFront).getPeakArea());//得到对应b ion位置的峰
            dMatch += Math.log10(lMSTwotrails.getPeakArea());///Math.log10(maxPeakArea);//(Math.log10(maxPeakArea)*iIonCurCount);//iIonCurCount;//得到对应b ion位置的峰


        }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
        {
            double minIonDistance = lMSTwotrails.getMzApex() - maxIonMaxMass/(1+Utils.thresholdPPM);

            dMatch *=Math.sqrt(1- minIonDistance*minIonDistance);//得到理论mz与实际质荷比距离

        }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
        {
            double minIonDistance = lMSTwotrails.getMzApex() - maxIonMaxMass/(1+Utils.thresholdPPM);

//            dMatch += Math.sqrt((1- minIonDistance*minIonDistance)* Math.log10(lMSTwotrails.getPeakArea())/(Math.log10(maxPeakArea)* iIonCurCount);//得到理论mz与实际质荷比距离，对应位置的峰面积iBionCurCount;//
//            dMatch += Math.sqrt(1- minIonDistance*minIonDistance)* Math.log10(lMSTwotrails.getPeakArea());
//            dMatch += Math.sqrt(1- minIonDistance*minIonDistance)* Math.log10(lMSTwotrails.getPeakArea())/iIonCurCount;
            dMatch += (1- minIonDistance)* Math.log10(lMSTwotrails.getPeakArea());///Math.log10(maxPeakArea);//(Math.log10(maxPeakArea)*iIonCurCount);//iIonCurCount;
        }
        return dMatch;
    }
    private static double getdMatchIonScore(List<MSTwoTrail> lMSTwotrails, double maxIonMaxMass, double dMatch, int iposFront, double maxPeakArea, int iIonCurCount) {
        if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
        {
            //                   dMatch += Math.log10(lMSTwotrails.get(ibposFront).getPeakArea())/(Math.log10(maxPeakArea)*iBionCurCount);//得到对应b ion位置的峰面积 iBionCurCount;//
//            dMatch += Math.log10(lMSTwotrails.get(iposFront).getPeakArea());//得到对应b ion位置的峰
            dMatch += Math.log10(lMSTwotrails.get(iposFront).getPeakArea());///Math.log10(maxPeakArea);//(Math.log10(maxPeakArea)*iIonCurCount);//iIonCurCount;//得到对应b ion位置的峰


        }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
        {
            double minIonDistance = lMSTwotrails.get(iposFront).getMzApex() - maxIonMaxMass/(1+Utils.thresholdPPM);

            dMatch *=Math.sqrt(1- minIonDistance*minIonDistance);//得到理论mz与实际质荷比距离

        }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
        {
            double minIonDistance = lMSTwotrails.get(iposFront).getMzApex() - maxIonMaxMass/(1+Utils.thresholdPPM);

//            dMatch += Math.sqrt((1- minIonDistance*minIonDistance)* Math.log10(lMSTwotrails.get(iposFront).getPeakArea())/(Math.log10(maxPeakArea)* iIonCurCount);//得到理论mz与实际质荷比距离，对应位置的峰面积iBionCurCount;//
//            dMatch += Math.sqrt(1- minIonDistance*minIonDistance)* Math.log10(lMSTwotrails.get(iposFront).getPeakArea());
//            dMatch += Math.sqrt(1- minIonDistance*minIonDistance)* Math.log10(lMSTwotrails.get(iposFront).getPeakArea())/iIonCurCount;
            dMatch += (1- minIonDistance)* Math.log10(lMSTwotrails.get(iposFront).getPeakArea());///Math.log10(maxPeakArea);//(Math.log10(maxPeakArea)*iIonCurCount);//iIonCurCount;
        }
        return dMatch;
    }
    public static  MSTwoTrail CompateTwoTrail(MSTwoTrail ms2Trail1,MSTwoTrail ms2Trail2,double dIonmass)
    {
        if(ms2Trail1==null && ms2Trail2==null)
        {
            return null;
        }else if(ms2Trail1==null)
        {
            return ms2Trail2;
        }else if(ms2Trail2==null)
        {
            return ms2Trail1;
        }else
        {
            if( Utils.IonMatchMode == Utils.IonMatchEnum.MaxAbundancesAreaInRange)
            {
                return ms2Trail1.getPeakArea()>ms2Trail2.getPeakArea()?ms2Trail1:ms2Trail2;

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.NearestInIon )
            {
                return (ms2Trail1.getMzApex()-dIonmass)<(ms2Trail2.getMzApex()-dIonmass)?ms2Trail1:ms2Trail2;

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.FirstInIonset )//第一个满足条件的值
            {
                return ms2Trail1;//返回当前时间点的第一个匹配Trail
            }
        }
        return null;
    }

    private static int getiIonCurCountWithCheckIsotope(List<MSTwoTrail> lMSTwotrails, int ifristpos, double IonMaxMass, List<Integer> ionpos) {
        int iposnext = ifristpos;
        int iCurCount = 0;//用于统计区间出现可匹配的离子数
        double dmaxPeakArea = 0.0;
        int iIonMaxPeakAreaPos = iposnext;
        double dMinDistance = 1000.0;
        int iIonNearestPos = iposnext;
        //若是在精度范围区间，开始计数
        for(; iposnext< lMSTwotrails.size() &&
                lMSTwotrails.get(iposnext).getMzApex()<= IonMaxMass; iposnext++)
        {
            // 找到每个离子的不同的匹配模式，可以用多态
            if( Utils.IonMatchMode == Utils.IonMatchEnum.MaxAbundancesAreaInRange)//最大丰度
            {
                if(lMSTwotrails.get(iposnext).getPeakArea()>dmaxPeakArea)
                {
                    dmaxPeakArea = lMSTwotrails.get(iposnext).getPeakArea();
                    iIonMaxPeakAreaPos = iposnext;
                }

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.NearestInIon )//离理论最近
            {
                if(lMSTwotrails.get(iposnext).getMzApex()-(IonMaxMass/(1+Utils.thresholdPPM))<dMinDistance)
                {
                    dMinDistance = lMSTwotrails.get(iposnext).getMzApex()-(IonMaxMass/(1+Utils.thresholdPPM));
                    iIonNearestPos = iposnext;
                }
            }

            iCurCount++;//判定离子匹配个数
        }

        if (iCurCount > 0) {//若有匹配个数，按不同的策略加入MS2trail的位置

            if( Utils.IonMatchMode == Utils.IonMatchEnum.MaxAbundancesAreaInRange)
            {

                ionpos.add(iIonMaxPeakAreaPos);

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.NearestInIon )
            {
                ionpos.add(iIonNearestPos);

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.FirstInIonset )//第一个满足条件的值
            {
                ionpos.add(ifristpos);

            }
        }
        return iCurCount;
    }


    private static int getiIonCurCountWithCheckIsotope(List<MSTwoTrail> lMSTwotrails, int ifristpos, double IonMaxMass,
//                                                       double dIonIsotopeMax, double dIonIsotopeMin,
                                                       double[][] arr_leftIsotope,
                                                       double[][] arr_rightIsotope, int iCharge,
                                                       List<Integer> ionpos) {
        int iposnext = ifristpos;
        int iCurCount = 0;//用于统计区间出现可匹配的离子数
        double dmaxPeakArea = 0.0;
        int iIonMaxPeakAreaPos = iposnext;
        double dMinDistance = 1000.0;
        int iIonNearestPos = iposnext;
        boolean bzIsIsotope = false;
        //若是在精度范围区间，开始计数
        for(; iposnext< lMSTwotrails.size() &&
                lMSTwotrails.get(iposnext).getMzApex()<= IonMaxMass; iposnext++)
        {
            // 找到每个离子的不同的匹配模式，可以用多态
            if( Utils.IonMatchMode == Utils.IonMatchEnum.MaxAbundancesAreaInRange)//最大丰度
            {
                if(lMSTwotrails.get(iposnext).getPeakArea()>dmaxPeakArea)
                {
                    dmaxPeakArea = lMSTwotrails.get(iposnext).getPeakArea();
                    iIonMaxPeakAreaPos = iposnext;
                }

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.NearestInIon )//离理论最近
            {
                if(Math.abs(lMSTwotrails.get(iposnext).getMzApex()-(IonMaxMass/(1+Utils.thresholdPPM)))<dMinDistance)
                {
                    dMinDistance = Math.abs(lMSTwotrails.get(iposnext).getMzApex()-(IonMaxMass/(1+Utils.thresholdPPM)));
                    iIonNearestPos = iposnext;
                }
            }

            iCurCount++;//判定离子匹配个数
        }

        if (iCurCount > 0) {//若有匹配个数，按不同的策略加入MS2trail的位置

            int imatchpos  = 0;
            if( Utils.IonMatchMode == Utils.IonMatchEnum.MaxAbundancesAreaInRange)
            {
                imatchpos = iIonMaxPeakAreaPos;
//                ionpos.add(iIonMaxPeakAreaPos);

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.NearestInIon )
            {
                imatchpos = iIonNearestPos;
//                ionpos.add(iIonNearestPos);

            }else if( Utils.IonMatchMode == Utils.IonMatchEnum.FirstInIonset )//第一个满足条件的值
            {
                imatchpos = ifristpos;
//                ionpos.add(ifristpos);

            }
            if (Utils.bzCheckMS2FragmentIsotope)
            {

//                for (int i = ifristpos;i>=0 &&
//                        lMSTwotrails.get(i).getMzApex()>= dIonIsotopeMin; i--)
//                {
//                    //mz小于1500，左侧峰大
//                    //mz大于1500，左侧峰大于右侧峰0.5即可
//                    if ( lMSTwotrails.get(i).getMzApex()<= dIonIsotopeMax &&
//                            (
//                             (lMSTwotrails.get(i).getMzApex()<1500 &&  lMSTwotrails.get(imatchpos).getPeakArea()< lMSTwotrails.get(i).getPeakArea())
//                             ||
//                             (lMSTwotrails.get(i).getMzApex()>=1500 &&  lMSTwotrails.get(imatchpos).getPeakArea()/2.0 < lMSTwotrails.get(i).getPeakArea())
//                        )
//                    )
//                    {
//                        //检测前面有isotope峰，该匹配峰无效
//                        bzIsIsotope = true;
//                        iCurCount = 0;
//                        break;
//                    }
//                }

//                boolean bzIsIsotopeBefore = false;
                int i = ifristpos;
                for (int ic = Utils.iMS2Charge;ic>0;ic--)//不管是哪个电荷的右侧isotope，都要排除
                {
                    for (;i>=0 &&
                            lMSTwotrails.get(i).getMzApex()>=  arr_leftIsotope[ic-1][0]; i--)
                    {
                        //mz小于1500，左侧峰大
                        //mz大于1500，左侧峰大于右侧峰0.5即可
                        if ( lMSTwotrails.get(i).getMzApex()<= arr_leftIsotope[ic-1][1] &&
                                (
                                        (lMSTwotrails.get(i).getMzApex()<1500 &&  lMSTwotrails.get(imatchpos).getPeakArea()< lMSTwotrails.get(i).getPeakArea())
                                                ||
                                                (lMSTwotrails.get(i).getMzApex()>=1500 &&  lMSTwotrails.get(imatchpos).getPeakArea()/2.0 < lMSTwotrails.get(i).getPeakArea())
                                )
                        )
                        {
                            //检测前面有isotope峰，该匹配峰无效
                            bzIsIsotope = true;
                            iCurCount = 0;
                            break;
                        }
                    }
                    if(bzIsIsotope) {
                        break;
                    }
                }

            }
            //判断是不是匹配错了离子峰，根据isotope的电荷
            if(bzIsIsotope==false && Utils.bzCheckMS2FragmentIsotopeCharge)
            {
                boolean bzHasIsotope = false;
                //电荷从大到小，则mz是从小到大
                for( int ic =Utils.iMS2Charge;ic>0;ic--)
                {
                    for(; iposnext< lMSTwotrails.size() &&
                            lMSTwotrails.get(iposnext).getMzApex()<= arr_rightIsotope[ic-1][1]; iposnext++)//[1]表示该电荷的isotope的mz最大ppm值、[0]表示最小值
                    {
                        if ( lMSTwotrails.get(iposnext).getMzApex()>= arr_rightIsotope[ic-1][0] &&
                                (
                                        //当前峰与右侧isotope峰比较，小于1500是当前峰的丰度要大于右侧峰的丰度
                                        (lMSTwotrails.get(iposnext).getMzApex()<1500 &&  lMSTwotrails.get(imatchpos).getPeakArea() > lMSTwotrails.get(iposnext).getPeakArea())
                                                ||
                                                //当前峰与右侧isotope峰比较，大于1500是当前峰的丰度要大于右侧峰的丰度的一半即可
                                                (lMSTwotrails.get(iposnext).getMzApex()>=1500 &&  lMSTwotrails.get(imatchpos).getPeakArea() > lMSTwotrails.get(iposnext).getPeakArea() /2.0)
                                )
                        )
                        {
                            bzHasIsotope = true;//检测了一个isotope
                            if( ic!=iCharge )//是其它电荷的isotope
                            {
                                bzIsIsotope = true;
                                iCurCount = 0;
                            }
                            break;
                        }
                    }
                    //只要检测一个即可
                    if (bzHasIsotope){
                        break;
                    }

                }
            }
            if ( bzIsIsotope == false )
            {
                ionpos.add(imatchpos);
            }
        }
        return iCurCount;
    }
/*
    public static MatchPeptdide PepMatch(Peptide pep, List<MSTwoTrail> lMSTwotrails,double maxPeakArea, Enums.ScoreMode scoreMode, BufferedWriter writer,int icount,double ms1rt,double ms1Mass,long ms1id,double ms1qs) {
        double dMatch = 1.0;
        if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances || Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
        {
            dMatch = 0.0;//add all abundance
        }
        int len = pep.b_ions.length;
//        int iCurrentPos = 0;
        bpos.clear();
        ypos.clear();

        bipos.clear();
        yipos.clear();


        bMatchDistance.clear();
        yMatchDistance.clear();
        int ibpos = 0;
        int iypos = 0;

        int ibnearpos = 0;
        int iynearpos = 0;

        int bMatchSize;
        int yMatchSize;

        for(int i = 0; i < len; ++i) {

            bMatchSize = bpos.size();
            yMatchSize = ypos.size();

            double minBIonDistance = Utils.thresholdPPM;
            double minYIonDistance = Utils.thresholdPPM;

            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( pep.b_ions[i]*(1 + Utils.thresholdPPM) > (lMSTwotrails.get(0).getMzApex()) )
                && ( pep.b_ions[i]*(1 - Utils.thresholdPPM) < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
                if(i>0)
                {
                    if(pep.b_ions[i]<pep.y_ions[len-i])//根据上一次的找的y ion 位置，更新b ion 需要找的空间
                        ibpos = binarySearch0(lMSTwotrails,ibpos,iypos,pep.b_ions[i]);
                    else
                        ibpos = binarySearch0(lMSTwotrails,iypos,lMSTwotrails.size(),pep.b_ions[i]);


                }else
                {
                    ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size(),pep.b_ions[i]);//全新查找
                }

//            ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size()-1,pep.b_ions[i]);
                ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibpos);
                if (ibnearpos < lMSTwotrails.size() - 1) {
                    //判断大于和小于b ion的相邻数，哪个离b ion 最近，且小于阈值

                    if (ibnearpos == 0) {
                        double dmassibnearposBack = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);
                        if (dmassibnearposBack < minBIonDistance * pep.b_ions[i]) {
                            minBIonDistance = Math.abs(dmassibnearposBack) / 1000000.0;
                            bpos.add(ibnearpos);
                            bMatchDistance.add(minBIonDistance);
                        }
                    } else {
                        double dmassibnearpos = Math.abs(lMSTwotrails.get(ibnearpos - 1).getMzApex() - pep.b_ions[i]);
                        double dmassibnearposBack = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);

                        if (dmassibnearpos < minBIonDistance * pep.b_ions[i]) {
                            if (dmassibnearposBack < dmassibnearpos) {
                                minBIonDistance = Math.abs(dmassibnearposBack) / 1000000.0;
                                bpos.add(ibnearpos);
                                bMatchDistance.add(minBIonDistance);
                            } else {
                                minBIonDistance = Math.abs(dmassibnearpos) / 1000000.0;
                                bpos.add(ibnearpos - 1);
                                bMatchDistance.add(minBIonDistance);

                            }
                        } else if (dmassibnearposBack < minBIonDistance * pep.b_ions[i]) {
                            minBIonDistance = Math.abs(dmassibnearposBack) / 1000000.0;
                            bpos.add(ibnearpos);
                            bMatchDistance.add(minBIonDistance);


                        }
                    }

                } else {
                    ibnearpos = lMSTwotrails.size() - 1;
                    double dmassibnearpos = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);
                    if (dmassibnearpos < minBIonDistance * pep.b_ions[i]) {
                        minBIonDistance = Math.abs(dmassibnearpos) / 1000000.0;
                        bpos.add(ibnearpos);
                        bMatchDistance.add(minBIonDistance);

                    }
                }
//            iCurrentPos = ibpos;
//            if(iCurrentPos<0) iCurrentPos=0;
                if (bpos.size() > 0)
                    ibpos = bpos.get(bpos.size() - 1);
                if (ibpos < 0)
                    ibpos = ibnearpos;
            }



            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( pep.y_ions[len-i-1]*(1 + Utils.thresholdPPM) > (lMSTwotrails.get(0).getMzApex()) )
                    && ( pep.y_ions[len-i-1]*(1 - Utils.thresholdPPM) < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
                //根据前面比较的结果，search 对应的y离子
                Boolean bzBYIonEqual = false;
                if (pep.b_ions[i] > pep.y_ions[len-i-1])//根据以前的找的b ion 位置，更新y ion 需要找的空间
                {
                    iypos = binarySearch0(lMSTwotrails,iypos,ibpos,pep.y_ions[len-i-1]);

                }else  if (pep.b_ions[i] < pep.y_ions[len-i-1])
                {
                    iypos = binarySearch0(lMSTwotrails,iypos,lMSTwotrails.size(),pep.y_ions[len-i-1]);
                }else//若b，y相等不需要找y ion
                {

                    minYIonDistance  = minBIonDistance;
                    ypos.add(ibnearpos);
                    bzBYIonEqual = true;
                    yMatchDistance.add(minYIonDistance);

                }
                if(!bzBYIonEqual)//不相等
                {
                    iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(),iypos);
                    if (iynearpos < lMSTwotrails.size() -1)
                    {
                        //判断大于和小于y ion的相邻数，哪个离Y ion 最近，且小于阈值

                        if(iynearpos==0)
                        {
                            double dmassiynearposBack = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[len-i-1]);
                            if (dmassiynearposBack < minYIonDistance * pep.y_ions[len-i-1])
                            {
                                minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                                ypos.add(iynearpos);
                                yMatchDistance.add(minYIonDistance);


                            }

                        }else
                        {
                            double dmassiynearpos = Math.abs(lMSTwotrails.get(iynearpos-1).getMzApex() - pep.y_ions[len-i-1]);
                            double dmassiynearposBack = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[len-i-1]);

                            if (dmassiynearpos < minYIonDistance * pep.y_ions[len-i-1])
                            {
                                if (dmassiynearposBack < dmassiynearpos)
                                {
                                    minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                                    ypos.add(iynearpos);
                                    yMatchDistance.add(minYIonDistance);
                                }else
                                {
                                    minYIonDistance  = Math.abs(dmassiynearpos)/1000000.0;
                                    ypos.add(iynearpos-1);
                                    yMatchDistance.add(minYIonDistance);

                                }
                            }else if (dmassiynearposBack < minYIonDistance * pep.y_ions[len-i-1])
                            {
                                minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                                ypos.add(iynearpos);
                                yMatchDistance.add(minYIonDistance);


                            }
                        }

                    }else
                    {
//                        System.out.println(lMSTwotrails.size()+"-lmstowtrlastelem:"+lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()
//                                +"-pepelement:"+pep.y_ions[len-i-1]+"-iypos:"+iypos+"--pos:"+iynearpos);
                        iynearpos = lMSTwotrails.size()-1;
                        double dmassiynearpos = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[len-i-1]);
                        if (dmassiynearpos < minYIonDistance * pep.y_ions[len-i-1])
                        {
                            minYIonDistance  = Math.abs(dmassiynearpos)/1000000.0;
                            ypos.add(iynearpos);
                            yMatchDistance.add(minYIonDistance);

                        }
                    }

                }
//            iCurrentPos = iypos;
//            if(iCurrentPos<0) iCurrentPos=0;
                if(ypos.size()>0)
                    iypos = ypos.get(ypos.size()-1);
                if(iypos < 0)
                    iypos = iynearpos;
            }


            if(bpos.size()>bMatchSize)
            {
                bipos.add(i);
                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
                {
                    dMatch += Math.log10(lMSTwotrails.get(ibpos).getPeakArea())/Math.log10(maxPeakArea);//得到对应b ion位置的峰面积
//                    dMatch += Math.log10(lMSTwotrails.get(ibpos).getPeaksSum());//得到对应b ion位置的峰

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
                {
                    dMatch*=(1-minBIonDistance);//得到理论mz与实际质荷比距离

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
                {
                    dMatch += (1-minBIonDistance)* Math.log10(lMSTwotrails.get(ibpos).getPeakArea())/Math.log10(maxPeakArea);//得到理论mz与实际质荷比距离，对应位置的峰面积

                }
            }
            if(ypos.size()>yMatchSize)
            {
                yipos.add(len-i-1);
                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
                {
                    dMatch += Math.log10(lMSTwotrails.get(iypos).getPeakArea())/Math.log10(maxPeakArea);//得到对应y ion位置的峰面积
//                   dMatch += Math.log10(lMSTwotrails.get(iypos).getPeaksSum());//得到对应b ion位置的峰

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
                {
                    dMatch*=(1-minYIonDistance);//得到理论mz与实际质荷比距离

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
                {
                    dMatch += (1-minYIonDistance)*Math.log10(lMSTwotrails.get(iypos).getPeakArea())/Math.log10(maxPeakArea);//得到理论mz与实际质荷比距离，对应位置的峰面积

                }
            }

//            for (MSTwoTrail msTwoTrail : lMSTwotrails) {
//
//                double[] mass = msTwoTrail.getMzs();
//
//                int iPosBegin = Arrays.binarySearch(mass,msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass);
//
//                iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);
//
//                for(double drt:mass) {
//                    if (Math.abs(drt-pep.b_ions[i]) < minBIonDistance * pep.b_ions[i])
//                    {
//                        minBIonDistance = Math.abs(drt-pep.b_ions[i])/1000000.0;//if has the less one
//                        bpos.add(i);
//                    }
//                    if (Math.abs(drt-pep.y_ions[i]) < minYIonDistance * pep.y_ions[i])
//                    {
//                        minYIonDistance = Math.abs(drt-pep.y_ions[i])/1000000.0;//if has the less one
//                        ypos.add(i);
//
//                    }
//
//                }
//                //if find yion or bions
//                if(bpos.contains(i)) dMatch+=(1-minBIonDistance);
//                if(ypos.contains(i)) dMatch+=(1-minYIonDistance);
//
//            }
        }
        MatchPeptdide matchPep = null;
        int bySum = bpos.size() + ypos.size();
        if(bySum>0)
        {
            if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance) {
                dMatch = dMatch * bySum / (lMSTwotrails.size());
            }
//            dMatch = dMatch*bySum/(2*len*lMSTwotrails.size());
            if(dMatch>0.0) {
//                try {
//                        writer.write(icount + "\t");

                    ////不打印
//                    String strDecoy="target";
//                    writer.write(pep.id + "\t" );
//                    if(pep.id.contains("Decoy"))
//                    {
//                        writer.write("decoy"+"\t");
//                        strDecoy = "decoy";
//                    }
//                    else
//                        writer.write("target"+"\t");
//                    writer.write( pep.composition + "\t" + pep.mass);
//
//                    writer.write( "\t"+ ms1Mass);
//                    writer.write("\t" + ms1rt);
//                    writer.write("\t" + lMSTwotrails.get(0).getRtApex());
//                    writer.write("\t" + bySum);
//                    writer.write("\t" + dMatch);
////
////
//                   writer.write("\t"+"B:" + bipos.toString());
////                    String strMatch = "[";
////                    for(double dm:bMatchDistance)
////                        strMatch+=dm*1000000+",";
////                    if (strMatch.length()>1)
////                        strMatch =  strMatch.substring(0,strMatch.length()-2)+"]";
////                    else
////                        strMatch += "]";
////
////
////
////                     writer.write("\t" + strMatch);
////
//                    writer.write("|Y:" + yipos.toString());
////                    writer.write("\t" + yMatchDistance.toString());
////
////                    strMatch = "[";
////                    for(double dm:yMatchDistance)
////                        strMatch+=dm*1000000+",";
////                    if (strMatch.length()>1)
////                        strMatch =  strMatch.substring(0,strMatch.length()-2)+"]";
////                    else
////                        strMatch += "]";
//
////                    String rt = "\t";
////                    for (int i : bpos) {
////                        rt += lMSTwotrails.get(i).getRtApex() + ",";
////                    }
////                    writer.write("\t" + rt);
////
////                    rt = "\t";
////                    for (int i : ypos) {
////                        rt += lMSTwotrails.get(i).getRtApex() + ",";
////                    }
////                    writer.write("\t" + rt + "\n");
//                    writer.write("\n");
////不打印




//
//                } catch (IOException e) {
//                    e.printStackTrace();
//                }
                String strDecoy="target";
                if(pep.id.contains("Decoy"))
                {
                    strDecoy = "decoy";
                }
                String mstwortapex =""+lMSTwotrails.get(0).getRtApex();
                String strb= "|B:" + bipos.stream().map(x->x+1).collect(Collectors.toList()).toString();
                String stry= "|Y:" + yipos.stream().map(x->x+1).collect(Collectors.toList()).toString();
                matchPep = new MatchPeptdide();
                matchPep.setValues(
                   ms1id, dMatch,pep.id,strDecoy,pep.composition,pep.mass,ms1Mass,ms1rt,ms1qs,
                                mstwortapex,"" + bySum,""+dMatch,strb
                                ,stry,bipos,yipos);

                    }

//                                                        System.out.println(icount.incrementAndGet() +"--------");
//                                                        System.out.println(peps.get(iPos).composition + peps.get(iPos).mass);
//                                                        System.out.println(msonetrail.getMass());
//                                                        System.out.println(score);

//                System.out.println(arrMSTwoRts[iRTPos]);

        }
        else
            dMatch = 0.0;//未匹配到

 //       return dMatch;
        return  matchPep ;

    }
*/
    public static double PepMatch2(Peptide pep, List<MSTwoTrail> lMSTwotrails, Enums.ScoreMode scoreMode) {
        double dMatch = 1.0;
        int len = pep.b_ions.length;
//        int iCurrentPos = 0;
        List<Integer> bpos = new ArrayList<>();
        List<Integer> ypos = new ArrayList<>();
        List<Double> bMatchDistance = new ArrayList<>();
        List<Double> yMatchDistance = new ArrayList<>();
        int ibpos = 0;
        int iypos = 0;

        int ibnearpos = 0;
        int iynearpos = 0;

        int bMatchSize;
        int yMatchSize;

        for(int i = 0; i < len; ++i) {

            bMatchSize = bpos.size();
            yMatchSize = ypos.size();

            double minBIonDistance = Utils.thresholdPPM;
            double minYIonDistance = Utils.thresholdPPM;

            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( pep.b_ions[i]*(1 + Utils.thresholdPPM) > (lMSTwotrails.get(0).getMzApex()) )
                    && ( pep.b_ions[i]*(1 - Utils.thresholdPPM) < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
                if(i>0)
                {
                    if(pep.b_ions[i]<pep.y_ions[len-i])//根据上一次的找的y ion 位置，更新b ion 需要找的空间
                    {
                        ibpos = binarySearch0(lMSTwotrails,ibpos,iypos,pep.b_ions[i]);
                    } else {
                        ibpos = binarySearch0(lMSTwotrails,iypos,lMSTwotrails.size(),pep.b_ions[i]);
                    }


                }else
                {
                    ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size(),pep.b_ions[i]);//全新查找
                }

//            ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size()-1,pep.b_ions[i]);
                ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(),ibpos);
                if (ibnearpos < lMSTwotrails.size() -1)
                {
                    //判断大于和小于b ion的相邻数，哪个离b ion 最近，且小于阈值

                    if(ibnearpos==0)
                    {
                        double dmassibnearposBack = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);
                        if (dmassibnearposBack < minBIonDistance * pep.b_ions[i])
                        {
                            minBIonDistance  = Math.abs(dmassibnearposBack)/1000000.0;
                            bpos.add(ibnearpos);
                            bMatchDistance.add(minBIonDistance);


                        }

                    }else
                    {
                        double dmassibnearpos = Math.abs(lMSTwotrails.get(ibnearpos-1).getMzApex() - pep.b_ions[i]);
                        double dmassibnearposBack = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);

                        if (dmassibnearpos < minBIonDistance * pep.b_ions[i])
                        {
                            if (dmassibnearposBack < dmassibnearpos)
                            {
                                minBIonDistance  = Math.abs(dmassibnearposBack)/1000000.0;
                                bpos.add(ibnearpos);
                                bMatchDistance.add(minBIonDistance);


                            }else
                            {
                                minBIonDistance  = Math.abs(dmassibnearpos)/1000000.0;
                                bpos.add(ibnearpos-1);
                                bMatchDistance.add(minBIonDistance);

                            }
                        }else if (dmassibnearposBack < minBIonDistance * pep.b_ions[i])
                        {
                            minBIonDistance  = Math.abs(dmassibnearposBack)/1000000.0;
                            bpos.add(ibnearpos);
                            bMatchDistance.add(minBIonDistance);


                        }
                    }

                }else
                {
                    ibnearpos = lMSTwotrails.size()-1;
                    double dmassibnearpos = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);
                    if (dmassibnearpos < minBIonDistance * pep.b_ions[i])
                    {
                        minBIonDistance  = Math.abs(dmassibnearpos)/1000000.0;
                        bpos.add(ibnearpos);
                        bMatchDistance.add(minBIonDistance);

                    }
                }
//            iCurrentPos = ibpos;
//            if(iCurrentPos<0) iCurrentPos=0;
                if(bpos.size()>0) {
                    ibpos = bpos.get(bpos.size()-1);
                }
                if(ibpos < 0) {
                    ibpos =ibnearpos;
                }
            }



            //跳过极端情况，少于lMSTwotrails第一个元素和大于lMSTwotrails的最后一个元素
            if (( pep.y_ions[len-i-1]*(1 + Utils.thresholdPPM) > (lMSTwotrails.get(0).getMzApex()) )
                    && ( pep.y_ions[len-i-1]*(1 - Utils.thresholdPPM) < (lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()) ))
            {
                //根据前面比较的结果，search 对应的y离子
                Boolean bzBYIonEqual = false;
                if (pep.b_ions[i] > pep.y_ions[len-i-1])//根据以前的找的b ion 位置，更新y ion 需要找的空间
                {
                    iypos = binarySearch0(lMSTwotrails,iypos,ibpos,pep.y_ions[len-i-1]);

                }else  if (pep.b_ions[i] < pep.y_ions[len-i-1])
                {
                    iypos = binarySearch0(lMSTwotrails,iypos,lMSTwotrails.size(),pep.y_ions[len-i-1]);
                }else//若b，y相等不需要找y ion
                {

                    minYIonDistance  = minBIonDistance;
                    ypos.add(ibnearpos);
                    bzBYIonEqual = true;
                    yMatchDistance.add(minYIonDistance);

                }
                if(!bzBYIonEqual)//不相等
                {
                    iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(),iypos);
                    if (iynearpos < lMSTwotrails.size() -1)
                    {
                        //判断大于和小于y ion的相邻数，哪个离Y ion 最近，且小于阈值

                        if(iynearpos==0)
                        {
                            double dmassiynearposBack = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[len-i-1]);
                            if (dmassiynearposBack < minYIonDistance * pep.y_ions[len-i-1])
                            {
                                minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                                ypos.add(iynearpos);
                                yMatchDistance.add(minYIonDistance);


                            }

                        }else
                        {
                            double dmassiynearpos = Math.abs(lMSTwotrails.get(iynearpos-1).getMzApex() - pep.y_ions[len-i-1]);
                            double dmassiynearposBack = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[len-i-1]);

                            if (dmassiynearpos < minYIonDistance * pep.y_ions[len-i-1])
                            {
                                if (dmassiynearposBack < dmassiynearpos)
                                {
                                    minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                                    ypos.add(iynearpos);
                                    yMatchDistance.add(minYIonDistance);
                                }else
                                {
                                    minYIonDistance  = Math.abs(dmassiynearpos)/1000000.0;
                                    ypos.add(iynearpos-1);
                                    yMatchDistance.add(minYIonDistance);

                                }
                            }else if (dmassiynearposBack < minYIonDistance * pep.y_ions[len-i-1])
                            {
                                minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                                ypos.add(iynearpos);
                                yMatchDistance.add(minYIonDistance);


                            }
                        }

                    }else
                    {
//                        System.out.println(lMSTwotrails.size()+"-lmstowtrlastelem:"+lMSTwotrails.get(lMSTwotrails.size()-1).getMzApex()
//                                +"-pepelement:"+pep.y_ions[len-i-1]+"-iypos:"+iypos+"--pos:"+iynearpos);
                        iynearpos = lMSTwotrails.size()-1;
                        double dmassiynearpos = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[len-i-1]);
                        if (dmassiynearpos < minYIonDistance * pep.y_ions[len-i-1])
                        {
                            minYIonDistance  = Math.abs(dmassiynearpos)/1000000.0;
                            ypos.add(iynearpos);
                            yMatchDistance.add(minYIonDistance);

                        }
                    }

                }
//            iCurrentPos = iypos;
//            if(iCurrentPos<0) iCurrentPos=0;
                if(ypos.size()>0) {
                    iypos = ypos.get(ypos.size()-1);
                }
                if(iypos < 0) {
                    iypos = iynearpos;
                }
            }


            if(bpos.size()>bMatchSize) {
                dMatch*=(1-minBIonDistance);
            }
            if(ypos.size()>yMatchSize) {
                dMatch*=(1-minYIonDistance);
            }

//            for (MSTwoTrail msTwoTrail : lMSTwotrails) {
//
//                double[] mass = msTwoTrail.getMzs();
//
//                int iPosBegin = Arrays.binarySearch(mass,msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass);
//
//                iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);
//
//                for(double drt:mass) {
//                    if (Math.abs(drt-pep.b_ions[i]) < minBIonDistance * pep.b_ions[i])
//                    {
//                        minBIonDistance = Math.abs(drt-pep.b_ions[i])/1000000.0;//if has the less one
//                        bpos.add(i);
//                    }
//                    if (Math.abs(drt-pep.y_ions[i]) < minYIonDistance * pep.y_ions[i])
//                    {
//                        minYIonDistance = Math.abs(drt-pep.y_ions[i])/1000000.0;//if has the less one
//                        ypos.add(i);
//
//                    }
//
//                }
//                //if find yion or bions
//                if(bpos.contains(i)) dMatch+=(1-minBIonDistance);
//                if(ypos.contains(i)) dMatch+=(1-minYIonDistance);
//
//            }
        }

        int bySum = bpos.size() + ypos.size();
        if(bySum>0) {
            dMatch = dMatch*bySum/(2*len);
        } else {
            dMatch = 0.0;//未匹配到
        }

        return dMatch;

    }


    public static double PepMatch1(Peptide pep, List<MSTwoTrail> lMSTwotrails, Enums.ScoreMode scoreMode) {
        double dMatch = 0.0;
        int len = pep.b_ions.length;
//        int iCurrentPos = 0;
        List<Integer> bpos = new ArrayList<>();
        List<Integer> ypos = new ArrayList<>();

        int ibpos = 0;
        int iypos = 0;

        int ibnearpos = 0;
        int iynearpos = 0;

        for(int i = 0; i < len; ++i) {

            double minBIonDistance = Utils.thresholdPPM;
            double minYIonDistance = Utils.thresholdPPM;

            ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size()-1,pep.b_ions[i]);
            ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(),ibpos);
            if (ibnearpos < lMSTwotrails.size() -1)
            {
                //判断大于和小于b ion的相邻数，哪个离b ion 最近，且小于阈值
                double dmassibnearpos = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);
                double dmassibnearposBack = Math.abs(lMSTwotrails.get(ibnearpos+1).getMzApex() - pep.b_ions[i]);

                if (dmassibnearpos < minBIonDistance * pep.b_ions[i])
                {
                    if (dmassibnearposBack < dmassibnearpos)
                    {
                        minBIonDistance  = Math.abs(dmassibnearposBack)/1000000.0;
                        bpos.add(ibnearpos+1);


                    }else
                    {
                        minBIonDistance  = Math.abs(dmassibnearpos)/1000000.0;
                        bpos.add(ibnearpos);

                    }
                }else if (dmassibnearposBack < minBIonDistance * pep.b_ions[i])
                {
                    minBIonDistance  = Math.abs(dmassibnearposBack)/1000000.0;
                    bpos.add(ibnearpos+1);


                }
            }else
            {
                double dmassibnearpos = Math.abs(lMSTwotrails.get(ibnearpos).getMzApex() - pep.b_ions[i]);
                if (dmassibnearpos < minBIonDistance * pep.b_ions[i])
                {
                    minBIonDistance  = Math.abs(dmassibnearpos)/1000000.0;
                    bpos.add(ibnearpos);

                }
            }
//            iCurrentPos = ibpos;
//            if(iCurrentPos<0) iCurrentPos=0;
            if(bpos.size()>0) {
                ibpos = bpos.get(bpos.size()-1);
            }
            if(ibpos < 0) {
                ibpos =0;
            }

            //根据前面比较的结果，search 对应的y离子
            Boolean bzBYIonEqual = false;
            if (pep.b_ions[i] > pep.y_ions[len-i-1])//根据以前的找的b ion 位置，更新y ion 需要找的空间
            {
                iypos = binarySearch0(lMSTwotrails,iypos,lMSTwotrails.size()-1,pep.y_ions[len-i-1]);

            }else  if (pep.b_ions[i] < pep.y_ions[len-i-1])
            {
                iypos = binarySearch0(lMSTwotrails,iypos,lMSTwotrails.size()-1,pep.y_ions[len-i-1]);
            }else//若b，y相等不需要找y ion
            {

                minYIonDistance  = minBIonDistance;
                ypos.add(ibnearpos);
                bzBYIonEqual = true;

            }
            if(!bzBYIonEqual)//不相等
            {
                iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(),iypos);
                if (iynearpos < lMSTwotrails.size() -1)
                {
                    //判断大于和小于y ion的相邻数，哪个离Y ion 最近，且小于阈值

                    double dmassiynearpos = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[i]);
                    double dmassiynearposBack = Math.abs(lMSTwotrails.get(iynearpos+1).getMzApex() - pep.y_ions[i]);

                    if (dmassiynearpos < minYIonDistance * pep.y_ions[i])
                    {
                        if (dmassiynearposBack < dmassiynearpos)
                        {
                            minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                            ypos.add(iynearpos+1);


                        }else
                        {
                            minYIonDistance  = Math.abs(dmassiynearpos)/1000000.0;
                            ypos.add(iynearpos);

                        }
                    }else if (dmassiynearposBack < minYIonDistance * pep.y_ions[i])
                    {
                        minYIonDistance  = Math.abs(dmassiynearposBack)/1000000.0;
                        ypos.add(iynearpos+1);


                    }
                }else
                {
                    double dmassiynearpos = Math.abs(lMSTwotrails.get(iynearpos).getMzApex() - pep.y_ions[i]);
                    if (dmassiynearpos < minYIonDistance * pep.y_ions[i])
                    {
                        minYIonDistance  = Math.abs(dmassiynearpos)/1000000.0;
                        ypos.add(iynearpos);

                    }
                }

            }
//            iCurrentPos = iypos;
//            if(iCurrentPos<0) iCurrentPos=0;
            if(ypos.size()>0) {
                iypos = ypos.get(ypos.size()-1);
            }
            if(iypos < 0) {
                iypos =0;
            }

            if(bpos.contains(i)) {
                dMatch+=(1-minBIonDistance);
            }
            if(ypos.contains(i)) {
                dMatch+=(1-minYIonDistance);
            }

//            for (MSTwoTrail msTwoTrail : lMSTwotrails) {
//
//                double[] mass = msTwoTrail.getMzs();
//
//                int iPosBegin = Arrays.binarySearch(mass,msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass);
//
//                iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);
//
//                for(double drt:mass) {
//                    if (Math.abs(drt-pep.b_ions[i]) < minBIonDistance * pep.b_ions[i])
//                    {
//                        minBIonDistance = Math.abs(drt-pep.b_ions[i])/1000000.0;//if has the less one
//                        bpos.add(i);
//                    }
//                    if (Math.abs(drt-pep.y_ions[i]) < minYIonDistance * pep.y_ions[i])
//                    {
//                        minYIonDistance = Math.abs(drt-pep.y_ions[i])/1000000.0;//if has the less one
//                        ypos.add(i);
//
//                    }
//
//                }
//                //if find yion or bions
//                if(bpos.contains(i)) dMatch+=(1-minBIonDistance);
//                if(ypos.contains(i)) dMatch+=(1-minYIonDistance);
//
//            }
        }

        int bySum = bpos.size() + ypos.size();
        if(bySum>0) {
            dMatch = dMatch/bySum;
        }

        return dMatch;

    }

    private static int binarySearch0(List<MSTwoTrail>  a, int fromIndex, int toIndex,
                                     double key) {
        int low = fromIndex;
        int high = toIndex - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            double midVal = a.get(mid).getMzApex();

            if (midVal < key) {
                low = mid + 1;  // Neither val is NaN, thisVal is smaller
            } else if (midVal > key) {
                high = mid - 1; // Neither val is NaN, thisVal is larger
            } else {
                long midBits = Double.doubleToLongBits(midVal);
                long keyBits = Double.doubleToLongBits(key);
                if (midBits == keyBits)     // Values are equal
                {
                    return mid;             // Key found
                } else if (midBits < keyBits) // (-0.0, 0.0) or (!NaN, NaN)
                {
                    low = mid + 1;
                } else                        // (0.0, -0.0) or (NaN, !NaN)
                {
                    high = mid - 1;
                }
            }
        }

        return -(low + 1);  // key not found.
    }
/*
    //     returns a list of timestamps
    public static DbMatch PepMatch(Peptide pep, IsolationWindow window, edu.uw.waterlooms.match.Enums.ScoreMode scoreMode) {
        pep.GenerateIons(); // generate y/b ions here
        ArrayList<Peak> rtSet = new ArrayList<>(); // set of tallest peaks in matched y/b ion trail
        int len = pep.b_ions.length;
        for(int i = 0; i < len; ++i) {
            // add to list; use integer as this is most memory-efficient
            ArrayList<Integer> trailY = window.FindTrails(pep.y_ions[i], 1);
            ArrayList<Integer> trailB = window.FindTrails(pep.b_ions[i], 1);
            if (trailY != null) {
                for (Integer j: trailY) {
                    // find place to insert
                    rtSet.add(Peak.getPeakFromIndex(j, window));
                }
            }
            // check if found
            if (trailB != null) {
                for (Integer j: trailB) {
                    rtSet.add(Peak.getPeakFromIndex(j, window));
                }
            }
        }

        // if no trails matches, return null
        if (rtSet.size() == 0) {
            return null;
        }

        // sort peaks by retention time to perform sliding window later
        Peak[] tsArray = rtSet.toArray(new Peak[rtSet.size()]);
        Arrays.sort(tsArray, Peak.compareByRt);

        ArrayList<MatchedInterval> scoreChunks = new ArrayList<>();
        double curTsStart = tsArray[0].rt; // denotes beginning of current chunk of ts
        int curIndex = 0; // keep track of current index
        while (curIndex != -1) {
            int nextStartIndex = -1; // index for start of next window
            double curTsEnd = curTsStart+ edu.uw.waterlooms.match.Config.tsThreshold; // timestamp of when current window ends
            double nextTsStart = curTsStart+ Config.tsThreshold/2; // timestamp of start of next window

            MatchedInterval interval = new MatchedInterval(0, new ArrayList<>());

            while(tsArray[curIndex].rt < curTsEnd){
                // updates nextStartIndex and curTsStart
                if (curIndex < tsArray.length-1 && tsArray[curIndex].rt < nextTsStart) {
                    if (tsArray[curIndex+1].rt >= nextTsStart) {
                        nextStartIndex = curIndex+1;
                        curTsStart=nextTsStart < tsArray[curIndex+1].rt? tsArray[curIndex+1].rt : nextTsStart;
                    }
                }

                // check if a similar peak already exists: selects first peak seen
                boolean exists = false;
                for (Peak p : interval.peaks) {
                    if (p.compareByMz(tsArray[curIndex]) == 0) {
                        exists = true;
                    }
                }

                if(!exists) {
                    if (scoreMode == Enums.ScoreMode.RELATIVE) {
//                        double maxRtRound = Utils.RoundToX(tsArray[curIndex].rt, Config.maxGap);
//                        double maxRtIntensity = IsolationWindowCollection.max.get(maxRtRound);
                        double maxRtRound = tsArray[curIndex].rt;
                        tsArray[curIndex].score = 100*tsArray[curIndex].score / window.max.get(maxRtRound);
                    }
                    // if scoreMode is absolute, the score is initialized as max intensity,
                    // logging both modes calculates correct score
                    interval.score+=EmpiricalScore(tsArray[curIndex].score);
                    interval.peaks.add(tsArray[curIndex]);
                }

                ++curIndex;
                if(curIndex == tsArray.length) {
                    break;
                }
            }

            scoreChunks.add(interval);
            curIndex = nextStartIndex;

        }

//        for (MatchedInterval cur: scoreChunks) {
//            System.out.println(cur.score);
//            System.out.println(cur.peaks.stream().map(Object::toString)
//                    .collect(Collectors.joining("\n")));
//
//            if(cur.peaks.get(0).rt >= 61 &&  cur.peaks.get(0).rt <= 61.5 ) {
//                int a = 0;
//            }
//        }

        // pick highest score chunk
        MatchedInterval ret = scoreChunks.get(0);
        for (MatchedInterval interval : scoreChunks) {
            if (interval.score > ret.score) {
                ret = interval;
            }
        }
        return new DbMatch(pep.id, pep.composition, ret);
    }
*/
}
