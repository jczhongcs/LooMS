package edu.uw.waterlooms.peptideMatch;

import org.apache.commons.math3.stat.regression.SimpleRegression;

import java.util.*;
import java.util.concurrent.atomic.AtomicReference;

public class MatchPeptdide implements Comparable<MatchPeptdide>{

    public long getlMsoneid() {
        return lMsoneid;
    }

    public long lMsoneid;

    public double getdMatch() {
        return dMatch;
    }

    public double dMatch;
    String strPepid;
    public String strDecoy;
    String strPepcomposition;
    double dPmass;
//    double dms1Mass;
//    double dms1rt;
//    double dms1QualityScore;

    public void setMsOneTrail(MSOneTrail msOneTrail) {
        this.msOneTrail = msOneTrail;
    }

    public MSOneTrail msOneTrail;
    public void setDms2rt(double dms2rt) {
        this.dms2rt = dms2rt;
    }

    double dms2rt;
    String strMS2rts;

    public void setStrMS1trailInMS2TrailCount(String strMS1trailInMS2TrailCount) {
        this.strMS1trailInMS2TrailCount = strMS1trailInMS2TrailCount;
    }

    String strMS1trailInMS2TrailCount;
//    String strdMatchs;
    String strbipos;
    String stryipos;
    public List<Integer> lBIonCovered;
    public List<Integer> lYIonCovered;
    public int iBRepeatCount;
    public int iYRepeatCount;

//     List<Integer> bpos ;
//     List<Integer> ypos;

    public void setMapBionMS2Trail(TreeMap<Integer, LinkedList<MSTwoTrail>> mapBionMS2Trail) {
        this.mapBionMS2Trail = mapBionMS2Trail;
    }

    public void setMapYionMS2Trail(TreeMap<Integer, LinkedList<MSTwoTrail>> mapYionMS2Trail) {
        this.mapYionMS2Trail = mapYionMS2Trail;
    }

    TreeMap<Integer, LinkedList<MSTwoTrail>> mapBionMS2Trail;
     TreeMap<Integer, LinkedList<MSTwoTrail>> mapYionMS2Trail;

    TreeMap<Integer, MSTwoTrail> mapCurRTBionMS2Trail;
    TreeMap<Integer, MSTwoTrail> mapCurRTYionMS2Trail;

//     List<Double> bMatchDistance ;
//     List<Double> yMatchDistance;

    public int iWindowSize;
    int iBConsecutiveCount;
    int iYConsecutiveCount;

    double dlogBionPeakAreaSum;
    double dlogYionPeakAreaSum;

    double dBionMassErrorSum;
    double dYionMassErrorSum;

    double dBionPeakHalfSum;
    double dYionPeakHalfSum;

    double dBionRetentionTimeErrorSum;
    double dYionRetentionTimeErrorSum;

    double dPeptideMutationRate;

    boolean bAdjustScore = false;

    double dPearsonCorrelationWithPredictMSMS;

    public void setPep(Peptide pep) {
        this.pep = pep;
    }

    public Peptide pep;

    public void setValues(long msoneid,double dmatch, String id, String strdecoy, String composition, double mass,
                         double ms1Mass, double ms1rt, double ms1qs,String s, String s1, String s2, String s3, String s4,
                          List<Integer> lBIon,
                          List<Integer> lYIon) {
        lMsoneid = msoneid;
        dMatch=dmatch;
        strPepid=id;
        strDecoy=strdecoy;
         strPepcomposition=composition;
         dPmass=mass;
//         dms1Mass=ms1Mass;
//         dms1rt=ms1rt;
//         dms1QualityScore = ms1qs;
//         strMS2rts=s;
//         strbySums=s1;
//         strdMatchs=s2;
         strbipos=s3;
         stryipos=s4;
         lBIonCovered = new ArrayList<>(lBIon);
         lYIonCovered = new ArrayList<>(lYIon);
         iBRepeatCount = 0;
         iYRepeatCount = 0;


    }

    public void combineResult(MatchPeptdide matchPep) {
        if(Utils.PeptideMatchMode==Utils.PeptideMatchEnum.COMBINEMODE)
        {
             dMatch+=matchPep.dMatch;

//             strMS2rts+=','+matchPep.strMS2rts;
//             strbySums+=','+matchPep.strbySums;
//             strdMatchs+=','+matchPep.strdMatchs;
             strbipos+=matchPep.strbipos;
             stryipos+=matchPep.stryipos;

             //calculate all repeat element for b y ion
            List<Integer> lBIonBac = new ArrayList<>(matchPep.lBIonCovered);//backup matchPep.lBIonCovered
            lBIonBac.retainAll(lBIonCovered);
            iBRepeatCount += lBIonBac.size();
            List<Integer> lYIonBac = new ArrayList<>(matchPep.lYIonCovered);//backup matchPep.lYIonCovered
            lYIonBac.retainAll(lYIonCovered);
            iYRepeatCount += lYIonBac.size();

             lBIonCovered.addAll(matchPep.lBIonCovered);
            lYIonCovered.addAll(matchPep.lYIonCovered);

        }
    }

    public void combineResultSetRepeatCoverCount(MatchPeptdide matchPep) {


//             strMS2rts+=','+matchPep.strMS2rts;
//             strbySums+=','+matchPep.strbySums;
//             strdMatchs+=','+matchPep.strdMatchs;
//            strbipos+=matchPep.strbipos;
//            stryipos+=matchPep.stryipos;

            //calculate all repeat element for b y ion
            List<Integer> lBIonBac = new ArrayList<>(matchPep.lBIonCovered);//backup matchPep.lBIonCovered
            lBIonBac.retainAll(lBIonCovered);
            iBRepeatCount += lBIonBac.size();
            List<Integer> lYIonBac = new ArrayList<>(matchPep.lYIonCovered);//backup matchPep.lYIonCovered
            lYIonBac.retainAll(lYIonCovered);
            iYRepeatCount += lYIonBac.size();

            lBIonCovered.addAll(matchPep.lBIonCovered);
            lYIonCovered.addAll(matchPep.lYIonCovered);


    }

    @Override
    public String toString() {
//        return "MatchPeptdide{" +
//                "lMsoneid=" + lMsoneid +
//                ", dMatch=" + dMatch + '\'' +
//                ", strPepid='" + strPepid + '\'' +
//                ", strDecoy='" + strDecoy + '\'' +
//                ", strPepcomposition='" + strPepcomposition + '\'' +
//                ", dPmass=" + dPmass +
//                ", dms1Mass=" + dms1Mass +
//                ", dms1rt=" + dms1rt +
//                ", strMS2rts='" + strMS2rts + '\'' +
//                ", strbySums='" + strbySums + '\'' +
//                ", strdMatchs='" + strdMatchs + '\'' +
//                ", strbipos='" + strbipos + '\'' +
//                ", stryipos='" + stryipos + '\'' +
//                '}';

        int bcount = (int) strbipos.chars().filter(ch -> ch == ',').count();
        int ycount = (int) stryipos.chars().filter(ch -> ch == ',').count();
        if(bcount==0)
        {
            if(strbipos.length()>5)
            {
                bcount = 1;
            }
        }else
        {
            bcount++;
        }
        if(ycount==0)
        {
            if(stryipos.length()>5)
            {
                ycount = 1;
            }
        }else
        {
            ycount++;
        }

        int iBCoveredCount =lBIonCovered.stream().distinct().toArray().length;
        int iYCoveredCount =lYIonCovered.stream().distinct().toArray().length;


        return "" +
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

//                '\t' + dms1Mass +
//                '\t' +dms1rt +
//                '\t' +dms1QualityScore +
                '\t' +dms2rt +
//                '\t' +strMS2rts +
//                '\t' +strbySums +
//                '\t' +strdMatchs +
                '\t' +strbipos +
                '\t' +stryipos +
                '\t' +dlogBionPeakAreaSum +
                '\t' +dlogYionPeakAreaSum +
                '\t' + dBionMassErrorSum +
                '\t' + dYionMassErrorSum +
                '\t' + dBionPeakHalfSum +
                '\t' + dYionPeakHalfSum +
                '\t' + (bcount) +//计算bion的个数
                '\t' + (ycount) +//计算yion的个数
                '\t' +strPepcomposition.length() +
                '\t' + (msOneTrail.getMass()-dPmass)/(msOneTrail.getMass()/1000000.0) +

//                '\t' + (dms1Mass-dPmass)/(dms1Mass/1000000.0) +
                '\t' + dPeptideMutationRate +
                '\t' + dBionRetentionTimeErrorSum +
                '\t' + dYionRetentionTimeErrorSum +
                '\t' + iBConsecutiveCount+
                '\t' + iYConsecutiveCount +
                '\t' + ((bAdjustScore)?(mapBionMS2Trail==null?"": "B:"+mapBionMS2Trail.keySet()+" Y:"+(mapYionMS2Trail==null?"": mapYionMS2Trail.keySet())):"")+
                '\t' + strMS1trailInMS2TrailCount +

//                '\t' + iBRepeatCount+
//                '\t' + iYRepeatCount +
//                '\t' + iBCoveredCount+
//                '\t' + iYCoveredCount +
//                '\t' + lBIonCovered.toString()+
//                '\t' + lYIonCovered.toString() +
                "";
    }

    @Override
    public int compareTo(MatchPeptdide o) {
        return Double.compare(dMatch,o.dMatch);
    }

    public void setPeakAreaMSErrorPeakHalfRank(double dlogBionPeakAreaSum, double dlogYionPeakAreaSum, double dBionMassErrorSum, double dYionMassErrorSum,
                                               double iBionPeakHalfSum, double iYionPeakHalfSum,
                                               double dBionRetentionTimeErrorSum, double dYionRetentionTimeErrorSum,
                                                double dPeptideMutationRate) {
        this.dlogBionPeakAreaSum = dlogBionPeakAreaSum;
        this.dlogYionPeakAreaSum = dlogYionPeakAreaSum;

        this.dBionMassErrorSum = dBionMassErrorSum;
        this.dYionMassErrorSum = dYionMassErrorSum;

        this.dBionPeakHalfSum = iBionPeakHalfSum;
        this.dYionPeakHalfSum = iYionPeakHalfSum;

        this.dBionRetentionTimeErrorSum = dBionRetentionTimeErrorSum;
        this.dYionRetentionTimeErrorSum = dYionRetentionTimeErrorSum;
        this.dPeptideMutationRate = dPeptideMutationRate;
    }

    public void setConsetivePeak(int iBConsective, int iYConsective) {
        iBConsecutiveCount = iBConsective;
        iYConsecutiveCount = iYConsective;
    }

    public void setRepeatCoverCount(MatchPeptdide matchPeptdideAll) {
        iBRepeatCount = matchPeptdideAll.iBRepeatCount;

        iYRepeatCount = matchPeptdideAll.iYRepeatCount;

        lBIonCovered = new ArrayList<>(matchPeptdideAll.lBIonCovered);
        lYIonCovered = new ArrayList<>(matchPeptdideAll.lYIonCovered);
    }

    public void calculateCombineScore() {
        int bcount = (int) strbipos.chars().filter(ch -> ch == ',').count();
        int ycount = (int) stryipos.chars().filter(ch -> ch == ',').count();
        if(bcount==0)
        {
            if(strbipos.length()>5)
            {
                bcount = 1;
            }
        }else
        {
            bcount++;
        }
        if(ycount==0)
        {
            if(stryipos.length()>5)
            {
                ycount = 1;
            }
        }else
        {
            ycount++;
        }
        dMatch = /*-0.5599952 + dlogBionPeakAreaSum *-0.0084543 + dlogYionPeakAreaSum * 0.0110678
                + dBionMassErrorSum *0.1406417 + dYionMassErrorSum *0.2640127
                + dBionPeakHalfSum *0.0600398 + dYionPeakHalfSum *0.2028273
                + bcount *0.0190089 + ycount * -0.0308959
                + strPepcomposition.length()  *-0.0051013 + (dms1Mass-dPmass)/(dms1Mass/1000000.0)  *0.0017562;
*/

        -0.3283891+
                dlogBionPeakAreaSum       *  0.0015489+
                dlogYionPeakAreaSum       *  0.0035140+
                dBionMassErrorSum         *  0.0815693+
                dYionMassErrorSum         *  0.1068236+
                dBionPeakHalfSum          *  0.0914548+
                dYionPeakHalfSum          *  0.0683649+
                bcount            * -0.1321523+
                ycount            * -0.0672249+
                strPepcomposition.length()             *  0.0060725+
                (msOneTrail.getMass()-dPmass)/(msOneTrail.getMass()/1000000.0)       *  0.0012327+
//                (dms1Mass-dPmass)/(dms1Mass/1000000.0)       *  0.0012327+
                dPeptideMutationRate      * -0.2018165+
                dBionRetentionTimeErrorSum*  0.0474480+
                dYionRetentionTimeErrorSum*  0.0866186+
                iBConsecutiveCount              *  0.1030586+
                iYConsecutiveCount              *  0.1367920;
    }

    public void addmapBionMS2Trail(Integer ibp, MSTwoTrail msTwoTrail) {
        if(mapBionMS2Trail==null)
        {
            mapBionMS2Trail = new TreeMap<>();
        }
        LinkedList<MSTwoTrail> listMS2Trail =  mapBionMS2Trail.get(ibp);
        if(listMS2Trail==null)
        {
            listMS2Trail = new LinkedList<MSTwoTrail>();
        }
        listMS2Trail.add(msTwoTrail);
        mapBionMS2Trail.put(ibp,listMS2Trail);
    }


    public void addmapYionMS2Trail(Integer ibp, MSTwoTrail msTwoTrail) {
        if(mapYionMS2Trail==null)
        {
            mapYionMS2Trail = new TreeMap<>();
        }
        LinkedList<MSTwoTrail> listMS2Trail =  mapYionMS2Trail.get(ibp);
        if(listMS2Trail==null)
        {
            listMS2Trail = new LinkedList<MSTwoTrail>();
        }
        listMS2Trail.add(msTwoTrail);
        mapYionMS2Trail.put(ibp,listMS2Trail);
    }


    public boolean adjustScore() {
        //若移除了Trail,NEED TO READJUST SCORE.
        int bcount = mapBionMS2Trail==null?0: mapBionMS2Trail.size();
        int ycount = mapYionMS2Trail==null?0: mapYionMS2Trail.size();
        int iBConsective = iBConsecutiveCount;
        int iYConsective = iYConsecutiveCount;
        boolean needAdjust = false;
        if (removeTrails(mapBionMS2Trail) )
        {

            needAdjust = true;
            iBConsective = 0;

            dlogBionPeakAreaSum = 0.0;
            dBionMassErrorSum = 0.0;
            dBionPeakHalfSum = 0.0;
            dBionRetentionTimeErrorSum = 0.0;



            Set<Integer> setBpos = mapBionMS2Trail.keySet();
            int iBeforeValue = -2;

            for(Integer iBP:setBpos)
            {
                if (iBP == iBeforeValue + 1) {
                    iBConsective++;
                }
                iBeforeValue = iBP;


                //用于输出中间值

                dlogBionPeakAreaSum += Math.log10(mapBionMS2Trail.get(iBP).get(0).getPeakArea());
                dBionMassErrorSum += 1.0 - Math.pow((((mapBionMS2Trail.get(iBP).get(0).getMzApex() - pep.b_ions[iBP]) / (mapBionMS2Trail.get(iBP).get(0).getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                dBionPeakHalfSum += Math.max(Math.log10(50.0 / mapBionMS2Trail.get(iBP).get(0).getPeakAreaHalfRank()), 0.0);
                dBionRetentionTimeErrorSum += 1 - Math.pow((mapBionMS2Trail.get(iBP).get(0).getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
            }

            iBConsecutiveCount = iBConsective;
        }
        if ( removeTrails(mapYionMS2Trail)) {

            needAdjust =true;
            iYConsective = 0;

            dlogYionPeakAreaSum = 0.0;
            dYionMassErrorSum = 0.0;
            dYionPeakHalfSum = 0.0;
            dYionRetentionTimeErrorSum = 0.0;

            Set<Integer> setYpos = mapYionMS2Trail.keySet();
            int iYeforeValue = -2;

            for(Integer iYP:setYpos)
            {
                if (iYP == iYeforeValue + 1) {
                    iYConsective++;
                }
                iYeforeValue = iYP;


                //用于输出中间值

                dlogYionPeakAreaSum += Math.log10(mapYionMS2Trail.get(iYP).get(0).getPeakArea());
                dYionMassErrorSum += 1.0 - Math.pow((((mapYionMS2Trail.get(iYP).get(0).getMzApex() - pep.y_ions[iYP]) / (mapYionMS2Trail.get(iYP).get(0).getMzApex() / 1000000.0)) / (Utils.thresholdPPM * 1000000.0)), 2);
                dYionPeakHalfSum += Math.max(Math.log10(50.0 / mapYionMS2Trail.get(iYP).get(0).getPeakAreaHalfRank()), 0.0);
                dYionRetentionTimeErrorSum += 1 - Math.pow((mapYionMS2Trail.get(iYP).get(0).getRtApex() - msOneTrail.getRt()) / (Utils.thresholdRT), 2);
            }

            iYConsecutiveCount = iYConsective;

        }

        if(needAdjust)
        {
            bAdjustScore = true;

            dMatch =

                    -0.3283891+
                            dlogBionPeakAreaSum       *  0.0015489+
                            dlogYionPeakAreaSum       *  0.0035140+
                            dBionMassErrorSum         *  0.0815693+
                            dYionMassErrorSum         *  0.1068236+
                            dBionPeakHalfSum          *  0.0914548+
                            dYionPeakHalfSum          *  0.0683649+
                            bcount            * -0.1321523+
                            ycount            * -0.0672249+
                            strPepcomposition.length()             *  0.0060725+
                            (msOneTrail.getMass()-dPmass)/(msOneTrail.getMass()/1000000.0)       *  0.0012327+
                            dPeptideMutationRate      * -0.2018165+
                            dBionRetentionTimeErrorSum*  0.0474480+
                            dYionRetentionTimeErrorSum*  0.0866186+
                            iBConsective              *  0.1030586+
                            iYConsective              *  0.1367920;
        }

//        SetMS2TrailSelected();
        return needAdjust;
    }

    public void setMS2TrailSelected() {
        if(mapBionMS2Trail!=null) {
            Set<Integer> setBpos = mapBionMS2Trail.keySet();

            for (Integer iBP : setBpos) {
                mapBionMS2Trail.get(iBP).get(0).setbSelected(1);
            }
        }
        if(mapYionMS2Trail!=null) {
            Set<Integer> setYpos = mapYionMS2Trail.keySet();

            for (Integer iBP : setYpos) {
                mapYionMS2Trail.get(iBP).get(0).setbSelected(1);
            }
        }
    }


    public void setMS2TrailSelectedSpeed() {
        if(mapBionMS2Trail!=null) {
            mapBionMS2Trail.forEach((k,v)->
            {
                v.get(0).setbSelected(1);
            }
        );
//
//            Set<Integer> setBpos = mapBionMS2Trail.keySet();
//
//            for (Integer iBP : setBpos) {
//                mapBionMS2Trail.get(iBP).get(0).setbSelected(1);
//            }
        }
        if(mapYionMS2Trail!=null) {
            mapYionMS2Trail.forEach((k,v)->
            {
                v.get(0).setbSelected(1);
            });
//            Set<Integer> setYpos = mapYionMS2Trail.keySet();
//
//            for (Integer iBP : setYpos) {
//                mapYionMS2Trail.get(iBP).get(0).setbSelected(1);
//            }
        }
    }

    public void setMS2TrailMAXSelectedSpeed() {
        if(mapBionMS2Trail!=null) {
            mapBionMS2Trail.forEach((k,v)->
                    {
                        v.get(0).setbSelected(3);
                    }
            );
//
//            Set<Integer> setBpos = mapBionMS2Trail.keySet();
//
//            for (Integer iBP : setBpos) {
//                mapBionMS2Trail.get(iBP).get(0).setbSelected(1);
//            }
        }
        if(mapYionMS2Trail!=null) {
            mapYionMS2Trail.forEach((k,v)->
            {
                v.get(0).setbSelected(3);
            });
//            Set<Integer> setYpos = mapYionMS2Trail.keySet();
//
//            for (Integer iBP : setYpos) {
//                mapYionMS2Trail.get(iBP).get(0).setbSelected(1);
//            }
        }
    }
    public Boolean needAdjustScore()
    {
        AtomicReference<Boolean> bzNeedAdjust = new AtomicReference<>(false);
        mapBionMS2Trail.forEach((x,y)-> {
            for(MSTwoTrail msTwoTrail : y)
            {
                if(msTwoTrail.isbSelected())
                {
                    bzNeedAdjust.set(true);
                    return ;
                }
            }
        });


        if(bzNeedAdjust.get())  return bzNeedAdjust.get();

        mapYionMS2Trail.forEach((x,y)-> {
            for(MSTwoTrail msTwoTrail : y)
            {
                if(msTwoTrail.isbSelected())
                {
                    bzNeedAdjust.set(true);
                    return ;
                }
            }
        });


       return bzNeedAdjust.get();
    }

    public boolean removeTrails(TreeMap<Integer, LinkedList<MSTwoTrail>> mapIonMS2Trail)
    {
        Boolean bzRemoveTrail = false;
        if (mapIonMS2Trail!=null) {
            Set<Integer> setDel = new TreeSet<>();
            Set<Integer> setBpos = mapIonMS2Trail.keySet();
            for (Integer iBP : setBpos) {
                LinkedList<MSTwoTrail> listMS2BionTrial = mapIonMS2Trail.get(iBP);
                for (int i = 0; i < listMS2BionTrial.size(); i++) {
                    if (listMS2BionTrial.get(i).isbSelected()) {
                        bzRemoveTrail = true;
                        listMS2BionTrial.remove(i);
                        i--;
                    }
                }
                if (listMS2BionTrial.size() == 0) {
                    setDel.add(iBP);
                }
            }
            for (Integer id : setDel) {
                bzRemoveTrail = true;

                mapIonMS2Trail.remove(id);
            }
        }
        return bzRemoveTrail;
    }



    public void calculatePearsonCorrelationWithPredictMSMS(double[][] arrPredictMSMS)
    {
        double dBionPearsonCor = 0.0;
        double dYionPearsonCor = 0.0;
        SimpleRegression regression = new SimpleRegression();

        if(mapBionMS2Trail!=null) {

            int iCharge = 1;
            for (Integer bi : mapBionMS2Trail.keySet()) {
                if (bi >= iCharge * strPepcomposition.length())//The next charge
                {
                    dBionPearsonCor += regression.getR();
                    iCharge++;
                    regression.clear();
                }

                regression.addData(mapBionMS2Trail.get(bi).get(0).getPeakArea(), arrPredictMSMS[iCharge-1][bi-iCharge * strPepcomposition.length()]);
            }

            dBionPearsonCor += regression.getR();//the lastTIME

        }


        if(mapYionMS2Trail!=null) {
            int iCharge = 1;

            regression.clear();

            for (Integer yi : mapYionMS2Trail.keySet()) {

                if (yi >= iCharge * strPepcomposition.length())//next charge
                {
                    dYionPearsonCor += regression.getR();
                    iCharge++;
                    regression.clear();
                }
                regression.addData(mapBionMS2Trail.get(yi).get(0).getPeakArea(), arrPredictMSMS[iCharge-1][yi-iCharge * strPepcomposition.length()]);

            }
            dYionPearsonCor += regression.getR();

        }
        dPearsonCorrelationWithPredictMSMS = dBionPearsonCor+dYionPearsonCor;

    }


}
