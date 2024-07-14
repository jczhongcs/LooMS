package edu.uw.waterlooms.peptideMatch;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class MatchPeptideWithMultiCharge extends  MatchPeptdide {
    public double dPrositRT;
    int[] arr_bSum;
    int[] arr_ySum;
    double[] arr_dlogBionPeakAreaSum;
    double[] arr_dlogYionPeakAreaSum;
    double[] arr_dBionMassErrorSum;
    double[] arr_dYionMassErrorSum;
    double[] arr_dBionPeakHalfSum;
    double[] arr_dYionPeakHalfSum;
    double[] arr_dBionRetentionTimeErrorSum;
    double[] arr_dYionRetentionTimeErrorSum;

    double[] arr_dBionPearsonCorrelationWithPredictMSMS;
    double[] arr_dYionPearsonCorrelationWithPredictMSMS;

    double dBion60PearsonCorrelationWithPredictMSMS;
    double dCosWithPredictMSMSMatchedPredict;
    double dCosWithPredictMSMSMatchedTop6AllPredict;
    double dCosWithPredictMSMSMatchedTop6Charge1Predict;
    double dCosWithPredictMSMSMatchedTop6BIonCharge1Predict;
    double dCosWithPredictMSMSMatchedTop6YIonCharge1Predict;

    //    double[] arr_dYion60PearsonCorrelationWithPredictMSMS;
    static SimilarityCalculate regressionAllStatic = new SimilarityCalculate(false);

    int[] arr_iBConsective;
    int[] arr_iYConsective;

    int ibAllMatchedWithloss = 0;
    int ibMatchNormalAndH2Oloss = 0;
    int ibMatchNormalAndNH3loss = 0;
    int iyAllMatchedWithloss = 0;
    int iyMatchNormalAndH2Oloss = 0;
    int iyMatchNormalAndNH3loss = 0;
    int ibyComplemetaryCount = 0;


    String strSharePeaksMaxInfo;

    String strSharePeaksAllInfo;


    public void setValues(double dmatch, String id, String strdecoy, String composition, double mass, MSOneTrail msOneTrail,
                          String mstwortapex, int[] bSum, int[] ySum, double dMatch1,
                          String strb, String stry, List<Integer> bipos, List<Integer> yipos) {
        lMsoneid = msOneTrail.getId();
        dMatch=dmatch;
        strPepid=id;
        strDecoy=strdecoy;
        strPepcomposition=composition;
        dPmass=mass;
        strbipos=strb;
        stryipos=stry;
//        lBIonCovered = new ArrayList<>(bipos);
//        lYIonCovered = new ArrayList<>(yipos);
        iBRepeatCount = 0;
        iYRepeatCount = 0;
        arr_bSum = bSum;
        arr_ySum = ySum;
    }

    public void setPeakAreaMSErrorPeakHalfRankWithMultiCharge(double[] dlogBionPeakAreaSum, double[] dlogYionPeakAreaSum, double[] dBionMassErrorSum,
                                                              double[] dYionMassErrorSum, double[] dBionPeakHalfSum, double[] dYionPeakHalfSum, 
                                                              double[] dBionRetentionTimeErrorSum, double[] dYionRetentionTimeErrorSum, double dMutationRate) {

       arr_dlogBionPeakAreaSum = dlogBionPeakAreaSum;
       arr_dlogYionPeakAreaSum = dlogYionPeakAreaSum;
       arr_dBionMassErrorSum = dBionMassErrorSum;
       arr_dYionMassErrorSum = dYionMassErrorSum;
       arr_dBionPeakHalfSum = dBionPeakHalfSum;
       arr_dYionPeakHalfSum = dYionPeakHalfSum;
       arr_dBionRetentionTimeErrorSum = dBionRetentionTimeErrorSum;
       arr_dYionRetentionTimeErrorSum = dYionRetentionTimeErrorSum;
       dPeptideMutationRate = dMutationRate;

    }

    public void setConsetivePeakWithMultiCharge(int[] iBConsective, int[] iYConsective) {
        arr_iBConsective = iBConsective;
        arr_iYConsective = iYConsective;
    }

    public void calculateCombineScore() {


        dMatch = /*-0.5599952 + dlogBionPeakAreaSum *-0.0084543 + dlogYionPeakAreaSum * 0.0110678
                + dBionMassErrorSum *0.1406417 + dYionMassErrorSum *0.2640127
                + dBionPeakHalfSum *0.0600398 + dYionPeakHalfSum *0.2028273
                + bcount *0.0190089 + ycount * -0.0308959
                + strPepcomposition.length()  *-0.0051013 + (dms1Mass-dPmass)/(dms1Mass/1000000.0)  *0.0017562;
*/

                -0.3283891+
                        arr_dlogBionPeakAreaSum[0]       *  0.0015489+
                        arr_dlogYionPeakAreaSum[0]       *  0.0035140+
                        arr_dBionMassErrorSum[0]         *  0.0815693+
                        arr_dYionMassErrorSum[0]         *  0.1068236+
                        arr_dBionPeakHalfSum[0]          *  0.0914548+
                        arr_dYionPeakHalfSum[0]          *  0.0683649+
                        arr_bSum[0]            * -0.1321523+
                        arr_ySum[0]            * -0.0672249+
                        strPepcomposition.length()             *  0.0060725+
                        (msOneTrail.getMass()-dPmass)/(msOneTrail.getMass()/1000000.0)       *  0.0012327+
//                (dms1Mass-dPmass)/(dms1Mass/1000000.0)       *  0.0012327+
                        dPeptideMutationRate      * -0.2018165+
                        arr_dBionRetentionTimeErrorSum[0]*  0.0474480+
                        arr_dYionRetentionTimeErrorSum[0]*  0.0866186+
                        arr_iBConsective[0]              *  0.1030586+
                        arr_iYConsective[0]              *  0.1367920;
    }

    @Override
    public String toString() {
/*        return "MatchPeptideWithMultiCharge{" +
                "arr_bSum=" + Arrays.toString(arr_bSum) +
                ", arr_ySum=" + Arrays.toString(arr_ySum) +
                ", arr_dlogBionPeakAreaSum=" + Arrays.toString(arr_dlogBionPeakAreaSum) +
                ", arr_dlogYionPeakAreaSum=" + Arrays.toString(arr_dlogYionPeakAreaSum) +
                ", arr_dBionMassErrorSum=" + Arrays.toString(arr_dBionMassErrorSum) +
                ", arr_dYionMassErrorSum=" + Arrays.toString(arr_dYionMassErrorSum) +
                ", arr_dBionPeakHalfSum=" + Arrays.toString(arr_dBionPeakHalfSum) +
                ", arr_dYionPeakHalfSum=" + Arrays.toString(arr_dYionPeakHalfSum) +
                ", arr_dBionRetentionTimeErrorSum=" + Arrays.toString(arr_dBionRetentionTimeErrorSum) +
                ", arr_dYionRetentionTimeErrorSum=" + Arrays.toString(arr_dYionRetentionTimeErrorSum) +
                ", arr_iBConsective=" + Arrays.toString(arr_iBConsective) +
                ", arr_iYConsective=" + Arrays.toString(arr_iYConsective) +
                '}';*/






        String str = "" +
                "" + lMsoneid +
                '\t' + dMatch +
                '\t' + strPepid +
                '\t' +strDecoy +
                '\t' + strPepcomposition +
                '\t' +dPmass +
                '\t' + msOneTrail.getMass() +
                '\t' +msOneTrail.getRt() +
                '\t' +msOneTrail.getQuality_score() +
                '\t' +msOneTrail.getMz() +
                '\t' +msOneTrail.getZ() +
                '\t' +msOneTrail.getdPeakAreaLocalRankInTimeSpan() +

                '\t' +dms2rt +
//                '\t' +strbipos +
//                '\t' +stryipos +
                '\t' + strPepcomposition.length() +
                '\t' + (msOneTrail.getMass() - dPmass) / (msOneTrail.getMass() / 1000000.0) +
                '\t' + dPeptideMutationRate ;

        String str1="";
        String strbi="";
        String stryi="";
        boolean bzcharge=true;
        boolean yzcharge=true;

        double dMS1WithMS2BTrailCoscinSimilarity = 0.0;
        double dMS1WithMS2YTrailCoscinSimilarity = 0.0;

        double MS2BTrailCoscinSimilarity = 0.0;
        double MS2YTrailCoscinSimilarity = 0.0;


        int imaxB1Num = 0;
        int imaxY1Num = 0;


        if(mapBionMS2Trail!=null)
        {

            int iCharge =  1;

            int iBeforeValue = -2;




            for(Integer bi:mapBionMS2Trail.keySet())
            {

                if (bi == iBeforeValue + 1) {//consective ion
                    MS2BTrailCoscinSimilarity += Utils.coscinSimilarity(mapBionMS2Trail.get(iBeforeValue).get(0).getRts(),mapBionMS2Trail.get(iBeforeValue).get(0).getInts(),
                            mapBionMS2Trail.get(bi).get(0).getRts(),mapBionMS2Trail.get(bi).get(0).getInts());

                }
                iBeforeValue = bi;
                if(bi<strPepcomposition.length()) imaxB1Num = bi+1;

                if(bi>=iCharge*strPepcomposition.length())//The next charge
                {
                    bzcharge = true;
                    iCharge++;
                    iBeforeValue = -2;

                }
                if(bzcharge )
                {

                    if (strbi.length()>1)
                        strbi =  strbi.substring(0,strbi.length()-1)+"]";

                    strbi=strbi+" "+(iCharge)+":[";
                    bzcharge = false;
                }

                strbi+=(bi+1-((iCharge-1)*strPepcomposition.length()))+",";


//                Arrays.sort(msOneTrail.getArrMS1Rts());

                dMS1WithMS2BTrailCoscinSimilarity += Utils.coscinSimilarity(msOneTrail.getArrMS1Rts(),msOneTrail.getArrMS1Ins(),
                        mapBionMS2Trail.get(bi).get(0).getRts(),mapBionMS2Trail.get(bi).get(0).getInts());
                /*int ims1 = 0;
                int ims2 = 0;
                double dotMultiply = 0.0;
                double normA = 0.0;
                double normB = 0.0;

                while (ims1<msOneTrail.getArrMS1Rts().length || ims2<mapBionMS2Trail.get(bi).get(0).getRts().length)
                {
                   for (;ims1==msOneTrail.getArrMS1Rts().length || mapBionMS2Trail.get(bi).get(0).getRts()[ims2]<msOneTrail.getArrMS1Rts()[ims1];ims2++)
                   {
                       dotMultiply += mapBionMS2Trail.get(bi).get(0).getRts()[ims2] *  0.0;
                       normA += Math.pow(mapBionMS2Trail.get(bi).get(0).getRts()[ims2], 2);
                       normB += 0.0;


                   }

                   for (;ims2==mapBionMS2Trail.get(bi).get(0).getRts().length
                           || (((ims1+1)<msOneTrail.getArrMS1Rts().length) && mapBionMS2Trail.get(bi).get(0).getRts()[ims2] > msOneTrail.getArrMS1Rts()[ims1 + 1]); ims1++)
                   {
                       dotMultiply += 0.0 * msOneTrail.getArrMS1Rts()[ims1];
                       normA += 0.0;
                       normB += Math.pow(msOneTrail.getArrMS1Rts()[ims1], 2);

                   }
                   dotMultiply += mapBionMS2Trail.get(bi).get(0).getRts()[ims2] * msOneTrail.getArrMS1Rts()[ims1];
                   normA += Math.pow(mapBionMS2Trail.get(bi).get(0).getRts()[ims2], 2);
                   normB += Math.pow(msOneTrail.getArrMS1Rts()[ims1], 2);
                   ims2++;
                   ims1++;
                }
                double consinSimilarity = dotMultiply / (Math.sqrt(normA) * Math.sqrt(normB));*/


                /*if(msOneTrail.getId()==24474) {
                    //calculate the shape of trail between ms1 and ms2
                    System.out.println("ms1rt:" + JSONWriter.valueToString(msOneTrail.getArrMS1Rts()));
                    System.out.println("ms1ins:" + JSONWriter.valueToString(msOneTrail.getArrMS1Ins()));

                    System.out.println("ms2rt:" + JSONWriter.valueToString(mapBionMS2Trail.get(bi).get(0).getRts()));
                    System.out.println("ms2ins:" + JSONWriter.valueToString(mapBionMS2Trail.get(bi).get(0).getInts()));
                    //calculate the shape of trail between neibor ms2
                }*/

            }
            if (strbi.length()>1)
                strbi =  strbi.substring(0,strbi.length()-1)+"]";
            else
                strbi += "]";
        }
        ibyComplemetaryCount = 0;
        String strByComplemetary = "";
        if(mapYionMS2Trail!=null)
        {
            int iCharge =  1;


            int iBeforeValue = -2;


            boolean first = true;


             for(Integer yi:mapYionMS2Trail.keySet())
            {

                if(first)
                {
                    if(yi<strPepcomposition.length())
                    {
                        imaxY1Num  = strPepcomposition.length()-yi;
                    }
                    first = false;
                }

                if (yi == iBeforeValue + 1) {//consective ion
                    MS2YTrailCoscinSimilarity += Utils.coscinSimilarity(mapYionMS2Trail.get(iBeforeValue).get(0).getRts(),mapYionMS2Trail.get(iBeforeValue).get(0).getInts()
                            ,mapYionMS2Trail.get(yi).get(0).getRts(),mapYionMS2Trail.get(yi).get(0).getInts());

                }
                iBeforeValue = yi;

                if(yi>=iCharge*strPepcomposition.length())//next charge
                {
                    yzcharge = true;
                    iCharge++;
                    iBeforeValue = -2;
                }
                if(yzcharge )
                {
                    if (stryi.length()>1)
                        stryi =  stryi.substring(0,stryi.length()-1)+"]";

                    stryi=stryi+" "+(iCharge)+":[";
                    yzcharge = false;
                }
//                stryi+=(((iCharge)*strPepcomposition.length())-yi+1)+",";
                stryi+=(((iCharge)*strPepcomposition.length())-yi)+",";

                dMS1WithMS2YTrailCoscinSimilarity += Utils.coscinSimilarity(msOneTrail.getArrMS1Rts(),msOneTrail.getArrMS1Ins(),mapYionMS2Trail.get(yi).get(0).getRts(),mapYionMS2Trail.get(yi).get(0).getInts());


                if(yi < strPepcomposition.length()) {//calculate  1 charge by complemetary
                    if (mapBionMS2Trail != null) {
                        for (Integer bi : mapBionMS2Trail.keySet()) {
                            if (bi >=  strPepcomposition.length()) {
                                break;
                            }
//                            if ((-yi + bi + 2) == 0) {//是否互补
                            if ((-yi + bi + 1) == 0) {//是否互补

                                    ibyComplemetaryCount++;
                                strByComplemetary +="b:"+(bi+1)+"_y:"+(strPepcomposition.length()-yi+1)+" ";
                            }
                        }
                    }
                }
            }
            if (stryi.length()>1)
                stryi =  stryi.substring(0,stryi.length()-1)+"]";
            else
                stryi += "]";
        }


        for(int i = 0;i<Utils.iMS2Charge;i++) {

            str1 = str1 +   '\t' + arr_dlogBionPeakAreaSum[i] +
                    '\t' + arr_dlogYionPeakAreaSum[i] +
                    '\t' + arr_dBionMassErrorSum[i] +
                    '\t' + arr_dYionMassErrorSum[i] +
                    '\t' + arr_dBionPeakHalfSum[i] +
                    '\t' + arr_dYionPeakHalfSum[i] +
                    '\t' + arr_bSum[i] +
                    '\t' + arr_ySum[i] +

                    '\t' + arr_dBionRetentionTimeErrorSum[i] +
                    '\t' + arr_dYionRetentionTimeErrorSum[i] +
                    '\t' + arr_iBConsective[i] +
                    '\t' + arr_iYConsective[i] +
                    '\t' + ((arr_dBionPearsonCorrelationWithPredictMSMS!=null)?arr_dBionPearsonCorrelationWithPredictMSMS[i]:"") +
                    '\t' + ((arr_dYionPearsonCorrelationWithPredictMSMS!=null)?arr_dYionPearsonCorrelationWithPredictMSMS[i]:"") ;
        }

        String strBMatched ="";//不输出细节
//        strBMatched = getIonMatchedWithPeakarea(mapBionMS2Trail);
//        strBMatched = getCurIonMatchedWithIntensity(mapCurRTBionMS2Trail);
        String strYMatched ="";//不输出细节
//        strYMatched = getIonMatchedWithPeakarea(mapYionMS2Trail);
//        strYMatched = getCurIonMatchedWithIntensity(mapCurRTYionMS2Trail);

//        if (mapBionMS2Trail!=null) {
//            for (Map.Entry<Integer, LinkedList<MSTwoTrail>> entry : mapBionMS2Trail.entrySet()) {
//                strYMatched = strYMatched + entry.getKey() + ",";
//                strYMatched = strYMatched + entry.getValue().get(0).getPeakArea() + ",";
//            }
//            strYMatched = strYMatched.substring(0,strYMatched.length()-1);
//        }


        str= str + '\t' +ibyComplemetaryCount+ /*" "+strByComplemetary+*/
                '\t' +dPearsonCorrelationWithPredictMSMS +
                '\t' +(dMS1WithMS2BTrailCoscinSimilarity+dMS1WithMS2YTrailCoscinSimilarity) +
                '\t' +(MS2BTrailCoscinSimilarity+MS2YTrailCoscinSimilarity) +
                '\t' +strbi +
                '\t' +stryi +
                str1+
                '\t' +strBMatched+
                '\t' +strYMatched+
                '\t' + ((bAdjustScore)?(mapBionMS2Trail==null?"": "B:"+mapBionMS2Trail.keySet()+" Y:"+(mapYionMS2Trail==null?"": mapYionMS2Trail.keySet())):"")+
                '\t' + strMS1trailInMS2TrailCount +
                '\t' + iWindowSize +
                '\t' + dBion60PearsonCorrelationWithPredictMSMS +
                '\t' + dCosWithPredictMSMSMatchedPredict +
                '\t' + dPrositRT +
                '\t' + ibAllMatchedWithloss +
                '\t' + ibMatchNormalAndH2Oloss +
                '\t' + ibMatchNormalAndNH3loss +
                '\t' + iyAllMatchedWithloss +
                '\t' + iyMatchNormalAndH2Oloss +
                '\t' + iyMatchNormalAndNH3loss +
                '\t' + (msOneTrail.isbAddPrecursor()?"1":"0") +
                '\t' + FastaFile.getFrequenceRateMultiple(strPepcomposition) +
                '\t' + imaxB1Num +
                '\t' + imaxY1Num +
                '\t' + dCosWithPredictMSMSMatchedTop6AllPredict +
                '\t' + dCosWithPredictMSMSMatchedTop6Charge1Predict +
                '\t' + dCosWithPredictMSMSMatchedTop6BIonCharge1Predict +
                '\t' + dCosWithPredictMSMSMatchedTop6YIonCharge1Predict +

                "";
        return str;
    }

    private String getCurIonMatchedWithIntensity(TreeMap<Integer, MSTwoTrail> mapCurRTionMS2Trail) {
        String  strMatched="";
        if (mapCurRTionMS2Trail !=null) {
            for (Map.Entry<Integer, MSTwoTrail> entry : mapCurRTionMS2Trail.entrySet()) {
                strMatched = strMatched + (entry.getKey() ) + ",";
                strMatched = strMatched + entry.getValue().getCurIntensity(dms2rt) + ",";
            }
            strMatched = strMatched.substring(0, strMatched.length()-1);
        }
        return strMatched;
    }

    private String getIonMatchedWithPeakarea(TreeMap<Integer, LinkedList<MSTwoTrail>> mapIonMS2Trail) {
        String  strMatched="";
        if (mapIonMS2Trail !=null) {
            for (Map.Entry<Integer, LinkedList<MSTwoTrail>> entry : mapIonMS2Trail.entrySet()) {
                strMatched = strMatched + (entry.getKey() +1) + ",";
                strMatched = strMatched + entry.getValue().get(0).getPeakArea() + ",";
            }
            strMatched = strMatched.substring(0, strMatched.length()-1);
        }
        return strMatched;
    }

    public boolean adjustScore() {
        return false;
    }

    public void calculatePearsonCorrelationWithPredictMSMS_JustAll(double[][] arrPredictMSMS)
    {
        regressionAllStatic.clear();


        arr_dBionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];
        arr_dYionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];

        if(mapBionMS2Trail!=null) {

            int iCharge = 1;
            for (Integer bi : mapBionMS2Trail.keySet()) {
                if (bi >= iCharge * strPepcomposition.length())//The next charge
                {
                    iCharge++;
                }
                if(bi-(iCharge-1) * strPepcomposition.length()< strPepcomposition.length()-1)
                {
                    regressionAllStatic.addData(mapBionMS2Trail.get(bi).get(0).getPeakArea(), arrPredictMSMS[bi-(iCharge-1) * strPepcomposition.length()][iCharge-1]);
                }
            }

        }

        if(mapYionMS2Trail!=null) {
            int iCharge = 1;


            for (Integer yi : mapYionMS2Trail.keySet()) {

                if (yi >= iCharge * strPepcomposition.length())//next charge
                {
                    iCharge++;
                }
                if(yi-(iCharge-1) * strPepcomposition.length() < strPepcomposition.length()-1)
                {
                    regressionAllStatic.addData(mapYionMS2Trail.get(yi).get(0).getPeakArea(), arrPredictMSMS[yi-(iCharge-1) * strPepcomposition.length()][2+iCharge-1]);
                }
            }

        }
        dPearsonCorrelationWithPredictMSMS = regressionAllStatic.getCosSimilarity();

    }


    public void calculatePearsonCorrelationWithPredictMSMS(double[][] arrPredictMSMS)
    {
        double dBionPearsonCor = 0.0;
        double dYionPearsonCor = 0.0;
        SimilarityCalculate regression = new SimilarityCalculate(false);
        SimilarityCalculate regressionAll = new SimilarityCalculate(false);

        arr_dBionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];
        arr_dYionPearsonCorrelationWithPredictMSMS = new double[Utils.iMS2Charge];

        if(mapBionMS2Trail!=null) {

            int iCharge = 1;
            for (Integer bi : mapBionMS2Trail.keySet()) {
                if (bi >= iCharge * strPepcomposition.length())//The next charge
                {
                    arr_dBionPearsonCorrelationWithPredictMSMS[iCharge-1] = regression.getCosSimilarity();

                    dBionPearsonCor += arr_dBionPearsonCorrelationWithPredictMSMS[iCharge-1];/*regression.getR();*/

                    iCharge++;
                    regression.clear();
                }
                if(bi-(iCharge-1) * strPepcomposition.length()< strPepcomposition.length()-1)
                {
                    regressionAll.addData(mapBionMS2Trail.get(bi).get(0).getPeakArea(), arrPredictMSMS[bi-(iCharge-1) * strPepcomposition.length()][iCharge-1]);
                    regression.addData(mapBionMS2Trail.get(bi).get(0).getPeakArea(), arrPredictMSMS[bi-(iCharge-1) * strPepcomposition.length()][iCharge-1]);
                }
            }
            arr_dBionPearsonCorrelationWithPredictMSMS[iCharge-1] = regression.getCosSimilarity();

            dBionPearsonCor +=arr_dBionPearsonCorrelationWithPredictMSMS[iCharge-1];/* regression.getR();//the lastTIME*/

        }


        if(mapYionMS2Trail!=null) {
            int iCharge = 1;

            regression.clear();

            for (Integer yi : mapYionMS2Trail.keySet()) {

                if (yi >= iCharge * strPepcomposition.length())//next charge
                {

                    arr_dYionPearsonCorrelationWithPredictMSMS[iCharge-1] = regression.getCosSimilarity();
                    dYionPearsonCor +=  arr_dYionPearsonCorrelationWithPredictMSMS[iCharge-1];/*regression.getR();*/
                    iCharge++;
                    regression.clear();
                }
                if(yi-(iCharge-1) * strPepcomposition.length() < strPepcomposition.length()-1)
                {
                    regressionAll.addData(mapYionMS2Trail.get(yi).get(0).getPeakArea(), arrPredictMSMS[yi-(iCharge-1) * strPepcomposition.length()][2+iCharge-1]);
                    regression.addData(mapYionMS2Trail.get(yi).get(0).getPeakArea(), arrPredictMSMS[yi-(iCharge-1) * strPepcomposition.length()][2+iCharge-1]);
                }

            }
            arr_dYionPearsonCorrelationWithPredictMSMS[iCharge-1] = regression.getCosSimilarity();

            dYionPearsonCor +=  arr_dYionPearsonCorrelationWithPredictMSMS[iCharge-1] ;/*regression.getR();*/

        }
        dPearsonCorrelationWithPredictMSMS = regressionAll.getCosSimilarity();
//        dPearsonCorrelationWithPredictMSMS = dBionPearsonCor+dYionPearsonCor;

    }


    public String getConciseString() {
//        msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsonemz\tmsonecharge\tBIon\tYIon"
        String str = "" +
                 lMsoneid +
                '\t' +dMatch+
                '\t' +strDecoy +
                '\t' + strPepcomposition +
                '\t' +dPmass +
                '\t' + msOneTrail.getMass() +
                '\t' +msOneTrail.getRt() +
                '\t' +msOneTrail.getQuality_score() +
                '\t' +msOneTrail.getMz() +
                '\t' +msOneTrail.getZ() +
                '\t' + getIonMatchedWithPeakarea(mapBionMS2Trail) +
                '\t' + getIonMatchedWithPeakarea(mapYionMS2Trail);
        return str;

    }


    public int getSharePeakWithOtherMatchPSMSpeed(MatchPeptideWithMultiCharge matchHighPSM) {

        AtomicInteger iSharePeaks= new AtomicInteger();
//        matchHighPSM.setMS2TrailSelected();//初始化，将以前的by离子匹配清为1.
        matchHighPSM.setMS2TrailSelectedSpeed();//初始化，将以前的by离子匹配清为1.
        TreeMap<Integer, LinkedList<MSTwoTrail>> maphighBionMS2Trail = matchHighPSM.mapBionMS2Trail;

        strSharePeaksMaxInfo = "";
        if(maphighBionMS2Trail!=null) {
            maphighBionMS2Trail.forEach((k,v)->
            {
                MSTwoTrail ms2Trail=v.get(0);
                if(ms2Trail.isbSelectedBY()==false)
                {
                    if (mapBionMS2Trail != null) {
                        mapBionMS2Trail.forEach((bi,v1)->
                        {
                            MSTwoTrail ms2TrailThis=v1.get(0);
                            if(ms2Trail==ms2TrailThis)
                            {
                                iSharePeaks.getAndIncrement();
                                strSharePeaksMaxInfo += "B" + bi + "_rt_" + ms2TrailThis.getRtApex() + "_mz_" + ms2TrailThis.getMzApex() + ",";
                                ms2Trail.setbSelected(2);

                            }
                        });
                    }
                    if (mapYionMS2Trail != null) {
                        mapYionMS2Trail.forEach((yi,v1)->
                        {
                            MSTwoTrail ms2TrailThis=v1.get(0);
                            if(ms2Trail==ms2TrailThis)
                            {
                                iSharePeaks.getAndIncrement();
                                strSharePeaksMaxInfo += "Y" + yi + "_rt_" + ms2TrailThis.getRtApex() + "_mz_" + ms2TrailThis.getMzApex() + ",";
                                ms2Trail.setbSelected(2);
                            }
                        });
                    }
                }
            });

//            for (Integer bhi : maphighBionMS2Trail.keySet()) {
//                if(maphighBionMS2Trail.get(bhi).get(0).isbSelectedBY()==false ) {
//                    if (mapBionMS2Trail != null) {
//                        for (Integer bi : mapBionMS2Trail.keySet()) {
//
//                            if (maphighBionMS2Trail.get(bhi).get(0) == mapBionMS2Trail.get(bi).get(0)) {
//                                iSharePeaks.getAndIncrement();
////                                strSharePeaksMaxInfo += "B" + bi + "_rt_" + mapBionMS2Trail.get(bi).get(0).getRtApex() + "_mz_" + mapBionMS2Trail.get(bi).get(0).getMzApex() + ",";
//                                maphighBionMS2Trail.get(bhi).get(0).setbSelected(2);
//                            }
//                            if(mapBionMS2Trail.get(bi).get(0).getMzApex()>maphighBionMS2Trail.get(bhi).get(0).getMzApex()) break;
//                        }
//                    }
//                    if (mapYionMS2Trail != null) {
//                        for (Integer yi : mapYionMS2Trail.keySet()) {
//
//
//                            if (maphighBionMS2Trail.get(bhi).get(0) == mapYionMS2Trail.get(yi).get(0)) {
//                                iSharePeaks.getAndIncrement();
////                                strSharePeaksMaxInfo += "Y" + yi + "_rt_" + mapYionMS2Trail.get(yi).get(0).getRtApex() + "_mz_" + mapYionMS2Trail.get(yi).get(0).getMzApex() + ",";
//                                maphighBionMS2Trail.get(bhi).get(0).setbSelected(2);
//
//                            }
//                            if(mapYionMS2Trail.get(yi).get(0).getMzApex()>maphighBionMS2Trail.get(bhi).get(0).getMzApex()) break;
//                        }
//                    }
//                }
//            }
        }


        TreeMap<Integer, LinkedList<MSTwoTrail>> maphighYionMS2Trail = matchHighPSM.mapYionMS2Trail;

        if(maphighYionMS2Trail!=null) {
            maphighYionMS2Trail.forEach((k,v)->
            {
                MSTwoTrail ms2Trail=v.get(0);
                if(ms2Trail.isbSelectedBY()==false)
                {
                    if (mapBionMS2Trail != null) {
                        mapBionMS2Trail.forEach((bi,v1)->
                        {
                            MSTwoTrail ms2TrailThis=v1.get(0);
                            if(ms2Trail==ms2TrailThis)
                            {
                                iSharePeaks.getAndIncrement();
                                strSharePeaksMaxInfo += "B" + bi + "_rt_" + ms2TrailThis.getRtApex() + "_mz_" + ms2TrailThis.getMzApex() + ",";
                                ms2Trail.setbSelected(2);
                            }
                        });
                    }
                    if (mapYionMS2Trail != null) {
                        mapYionMS2Trail.forEach((yi,v1)->
                        {
                            MSTwoTrail ms2TrailThis=v1.get(0);
                            if(ms2Trail==ms2TrailThis)
                            {
                                iSharePeaks.getAndIncrement();
                                strSharePeaksMaxInfo += "Y" + yi + "_rt_" + ms2TrailThis.getRtApex() + "_mz_" + ms2TrailThis.getMzApex() + ",";
                                ms2Trail.setbSelected(2);
                            }
                        });
                    }
                }
            });

//            for (Integer hi : maphighYionMS2Trail.keySet()) {
//                if(maphighYionMS2Trail.get(hi).get(0).isbSelectedBY()==false ) {
//                    if (mapBionMS2Trail != null) {
//                        for (Integer bi : mapBionMS2Trail.keySet()) {
//
//                            if (maphighYionMS2Trail.get(hi).get(0) == mapBionMS2Trail.get(bi).get(0)) {
//                                iSharePeaks.getAndIncrement();
////                                strSharePeaksMaxInfo += "B" + bi + "_rt_" + mapBionMS2Trail.get(bi).get(0).getRtApex() + "_mz_" + mapBionMS2Trail.get(bi).get(0).getMzApex() + ",";
//                                maphighYionMS2Trail.get(hi).get(0).setbSelected(2);
//
//                            }
//                            if(mapBionMS2Trail.get(bi).get(0).getMzApex()>maphighYionMS2Trail.get(hi).get(0).getMzApex()) break;
//
//                        }
//                    }
//                    if (mapYionMS2Trail != null) {
//                        for (Integer yi : mapYionMS2Trail.keySet()) {
//
//                            if (maphighYionMS2Trail.get(hi).get(0) == mapYionMS2Trail.get(yi).get(0)) {
//                                iSharePeaks.getAndIncrement();
////                                strSharePeaksMaxInfo += "Y" + yi + "_rt_" + mapYionMS2Trail.get(yi).get(0).getRtApex() + "_mz_" + mapYionMS2Trail.get(yi).get(0).getMzApex() + ",";
//                                maphighYionMS2Trail.get(hi).get(0).setbSelected(2);
//
//                            }
//                            if(mapYionMS2Trail.get(yi).get(0).getMzApex()>maphighYionMS2Trail.get(hi).get(0).getMzApex()) break;
//
//                        }
//                    }
//                }
//            }
        }



        return iSharePeaks.get();
    }

    public int getSharePeakWithOtherMatchPSM(MatchPeptideWithMultiCharge matchHighPSM) {

        int iSharePeaks= 0;
//        matchHighPSM.setMS2TrailSelected();//初始化，将以前的by离子匹配清为1.
        matchHighPSM.setMS2TrailSelectedSpeed();//初始化，将以前的by离子匹配清为1.
        TreeMap<Integer, LinkedList<MSTwoTrail>> maphighBionMS2Trail = matchHighPSM.mapBionMS2Trail;

        strSharePeaksMaxInfo = "";
        if(maphighBionMS2Trail!=null) {
            for (Integer bhi : maphighBionMS2Trail.keySet()) {
                if(maphighBionMS2Trail.get(bhi).get(0).isbSelectedBY()==false ) {
                    if (mapBionMS2Trail != null) {
                        for (Integer bi : mapBionMS2Trail.keySet()) {

                            if (maphighBionMS2Trail.get(bhi).get(0) == mapBionMS2Trail.get(bi).get(0)) {
                                iSharePeaks++;
                                strSharePeaksMaxInfo += "B" + bi + "_rt_" + mapBionMS2Trail.get(bi).get(0).getRtApex() + "_mz_" + mapBionMS2Trail.get(bi).get(0).getMzApex() + ",";
                                maphighBionMS2Trail.get(bhi).get(0).setbSelected(2);
                            }
                        }
                    }
                    if (mapYionMS2Trail != null) {
                        for (Integer yi : mapYionMS2Trail.keySet()) {


                            if (maphighBionMS2Trail.get(bhi).get(0) == mapYionMS2Trail.get(yi).get(0)) {
                                iSharePeaks++;
                                strSharePeaksMaxInfo += "Y" + yi + "_rt_" + mapYionMS2Trail.get(yi).get(0).getRtApex() + "_mz_" + mapYionMS2Trail.get(yi).get(0).getMzApex() + ",";
                                maphighBionMS2Trail.get(bhi).get(0).setbSelected(2);

                            }
                        }
                    }
                }
            }
        }


        TreeMap<Integer, LinkedList<MSTwoTrail>> maphighYionMS2Trail = matchHighPSM.mapYionMS2Trail;

        if(maphighYionMS2Trail!=null) {

            for (Integer hi : maphighYionMS2Trail.keySet()) {
                if(maphighYionMS2Trail.get(hi).get(0).isbSelectedBY()==false ) {
                    if (mapBionMS2Trail != null) {
                        for (Integer bi : mapBionMS2Trail.keySet()) {

                            if (maphighYionMS2Trail.get(hi).get(0) == mapBionMS2Trail.get(bi).get(0)) {
                                iSharePeaks++;
                                strSharePeaksMaxInfo += "B" + bi + "_rt_" + mapBionMS2Trail.get(bi).get(0).getRtApex() + "_mz_" + mapBionMS2Trail.get(bi).get(0).getMzApex() + ",";
                                maphighYionMS2Trail.get(hi).get(0).setbSelected(2);

                            }
                        }
                    }
                    if (mapYionMS2Trail != null) {
                        for (Integer yi : mapYionMS2Trail.keySet()) {

                            if (maphighYionMS2Trail.get(hi).get(0) == mapYionMS2Trail.get(yi).get(0)) {
                                iSharePeaks++;
                                strSharePeaksMaxInfo += "Y" + yi + "_rt_" + mapYionMS2Trail.get(yi).get(0).getRtApex() + "_mz_" + mapYionMS2Trail.get(yi).get(0).getMzApex() + ",";
                                maphighYionMS2Trail.get(hi).get(0).setbSelected(2);

                            }
                        }
                    }
                }
            }
        }



        return iSharePeaks;
    }

    public int getSharePeakALLSpeed() {

        strSharePeaksAllInfo ="";
        AtomicInteger iSharePeaksALL = new AtomicInteger();
        if(mapBionMS2Trail!=null) {
            strSharePeaksAllInfo+="B:";

            mapBionMS2Trail.forEach((iBP,v)->
            {
                MSTwoTrail ms2Trail=v.get(0);

                if(ms2Trail.isbSelectedALL()) {
                    iSharePeaksALL.getAndIncrement();
                    strSharePeaksAllInfo+=iBP+"_rt_"+ms2Trail.getRtApex()+"_mz_"+ms2Trail.getMzApex()+",";
                }
            });
//            Set<Integer> setBpos = mapBionMS2Trail.keySet();
//
//            strSharePeaksAllInfo+="B:";
//            for (Integer iBP : setBpos) {
//                if(mapBionMS2Trail.get(iBP).get(0).isbSelected()) {
//                    iSharePeaksALL.getAndIncrement();
//                    strSharePeaksAllInfo+=iBP+"_rt_"+mapBionMS2Trail.get(iBP).get(0).getRtApex()+"_mz_"+mapBionMS2Trail.get(iBP).get(0).getMzApex()+",";
//                }
//            }
        }
        if(mapYionMS2Trail!=null) {
            strSharePeaksAllInfo+="Y:";
            mapYionMS2Trail.forEach((iYP,v)->
            {
                MSTwoTrail ms2Trail=v.get(0);

                if(ms2Trail.isbSelectedALL()) {
                    iSharePeaksALL.getAndIncrement();
                    strSharePeaksAllInfo+=iYP+"_rt_"+ms2Trail.getRtApex()+"_mz_"+ms2Trail.getMzApex()+",";
                }
            });

//            Set<Integer> setYpos = mapYionMS2Trail.keySet();
//            strSharePeaksAllInfo+="Y:";
//
//            for (Integer iBP : setYpos) {
//                if(mapYionMS2Trail.get(iBP).get(0).isbSelected()) {
//                    iSharePeaksALL.getAndIncrement();
//                    strSharePeaksAllInfo+=iBP+"_rt_"+mapYionMS2Trail.get(iBP).get(0).getRtApex()+"_mz_"+mapYionMS2Trail.get(iBP).get(0).getMzApex()+",";
//
//                }
//            }
        }
        return iSharePeaksALL.get();
    }

    public int getSharePeakALL() {

        strSharePeaksAllInfo ="";
        int iSharePeaksALL = 0;
        if(mapBionMS2Trail!=null) {
            Set<Integer> setBpos = mapBionMS2Trail.keySet();

            strSharePeaksAllInfo+="B:";
            for (Integer iBP : setBpos) {
                if(mapBionMS2Trail.get(iBP).get(0).isbSelected()) {
                    iSharePeaksALL++ ;
                    strSharePeaksAllInfo+=iBP+"_rt_"+mapBionMS2Trail.get(iBP).get(0).getRtApex()+"_mz_"+mapBionMS2Trail.get(iBP).get(0).getMzApex()+",";
                }
            }
        }
        if(mapYionMS2Trail!=null) {
            Set<Integer> setYpos = mapYionMS2Trail.keySet();
            strSharePeaksAllInfo+="Y:";

            for (Integer iBP : setYpos) {
                if(mapYionMS2Trail.get(iBP).get(0).isbSelected()) {
                    iSharePeaksALL++;
                    strSharePeaksAllInfo+=iBP+"_rt_"+mapYionMS2Trail.get(iBP).get(0).getRtApex()+"_mz_"+mapYionMS2Trail.get(iBP).get(0).getMzApex()+",";

                }
            }
        }
        return iSharePeaksALL;
    }

    public int getSharePeakWithOtherMatchWithNumberPSMSpeed() {
        strSharePeaksMaxInfo ="";
        AtomicInteger iSharePeaksMax = new AtomicInteger();
        if(mapBionMS2Trail!=null) {
            strSharePeaksMaxInfo+="B:";

            mapBionMS2Trail.forEach((iBP,v)->
            {
                MSTwoTrail ms2Trail=v.get(0);

                if(ms2Trail.isbSelectedMax()) {
                    iSharePeaksMax.getAndIncrement();
                    strSharePeaksMaxInfo+=iBP+"_rt_"+ms2Trail.getRtApex()+"_mz_"+ms2Trail.getMzApex()+",";
//                    ms2Trail.setbSelected(2);//不重复选

                }
            });

        }
        if(mapYionMS2Trail!=null) {
            strSharePeaksMaxInfo+="Y:";
            mapYionMS2Trail.forEach((iYP,v)->
            {
                MSTwoTrail ms2Trail=v.get(0);

                if(ms2Trail.isbSelectedMax()) {
                    iSharePeaksMax.getAndIncrement();
                    strSharePeaksMaxInfo+=iYP+"_rt_"+ms2Trail.getRtApex()+"_mz_"+ms2Trail.getMzApex()+",";
                }
            });


        }
        return iSharePeaksMax.get();
    }


}
