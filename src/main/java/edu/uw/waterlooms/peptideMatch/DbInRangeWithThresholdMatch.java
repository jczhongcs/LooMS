package edu.uw.waterlooms.peptideMatch;

import java.io.BufferedWriter;
import java.util.List;
import java.util.stream.Collectors;

public class DbInRangeWithThresholdMatch extends DbBaseMatch{



    @Override
    public  MatchPeptdide PepMatch(Peptide pep, List<MSTwoTrail> lMSTwotrails, double maxPeakArea,
                                   Enums.ScoreMode scoreMode, BufferedWriter writer, int icount,
                                   double ms1rt, double ms1Mass, long ms1id,double ms1qs) {

//        public DbBaseMatch PepMatch(Peptide pep, List<MSTwoTrail> lMSTwotrails,double maxPeakArea) {
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


        int ibposFront = 0;
        int ibposBack = 0;
        int iyposFront = 0;
        int iyposBack = 0;

        int ibnearpos = 0;
        int iynearpos = 0;

        int bMatchSize;
        int yMatchSize;

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
                if(i>0)
                {
                    if(bIonMinusMass<pep.y_ions[len-i])//根据上一次的找的y ion 位置，更新b ion 需要找的空间
                    {
                        ibposFront = binarySearch0(lMSTwotrails,ibposFront,iyposFront,bIonMinusMass);
                    }
                    else
                    {
                        ibposFront = binarySearch0(lMSTwotrails,iyposFront,lMSTwotrails.size(),bIonMinusMass);
                    }


                }else
                {
                    ibposFront = binarySearch0(lMSTwotrails,ibposFront,lMSTwotrails.size(),bIonMinusMass);//全新查找
                }

//            ibpos = binarySearch0(lMSTwotrails,ibpos,lMSTwotrails.size()-1,pep.b_ions[i]);
                ibnearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposFront);
                ibposBack = binarySearch0(lMSTwotrails,ibnearpos,lMSTwotrails.size(),bIonMaxMass);
                if (ibposBack < 0)//在区间内
                {
                    ibposBack = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), ibposBack);
                    if (ibposBack>0)
                        ibposBack--;
                }

                 iBionCurCount = ibposBack - ibnearpos;
                if( iBionCurCount>0 )
                {
                    bpos.add(ibnearpos);
                }
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
                Boolean bzBYIonEqual = false;
                if (bIonMinusMass > yIonMinusMass)//根据以前的找的b ion 位置，更新y ion 需要找的空间
                {
                    iyposFront = binarySearch0(lMSTwotrails,iyposFront,ibposFront,yIonMinusMass);

                }else  if (bIonMinusMass < yIonMinusMass)
                {
                    iyposFront = binarySearch0(lMSTwotrails,iyposFront,lMSTwotrails.size(),yIonMinusMass);
                }else//若b，y相等不需要找y ion
                {

                    minYIonDistance  = minBIonDistance;
                    ypos.add(ibnearpos);
                    bzBYIonEqual = true;
//                    yMatchDistance.add(minYIonDistance);

                }
                if(!bzBYIonEqual)//不相等
                {
                    iynearpos = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposFront);
                    iyposBack = binarySearch0(lMSTwotrails, iynearpos, lMSTwotrails.size(), yIonMaxMass);
                    if (iyposBack < 0)//在区间内
                    {
                        iyposBack = Utils.getiPosBeginOfNearest(lMSTwotrails.size(), iyposBack);
                        if (iyposBack > 0)
                            iyposBack--;
                    }

                    iYionCurCount = iyposBack - iynearpos;
                    if (iYionCurCount > 0) {
                        ypos.add(iynearpos);
                    }

                }
                if(ypos.size()>0)
                    iyposFront = ypos.get(ypos.size()-1);
                if(iyposFront < 0)
                    iyposFront = iynearpos;
            }


            if(bpos.size()>bMatchSize)
            {
                bipos.add(i);
                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
                {
                    dMatch += Math.log10(lMSTwotrails.get(ibposFront).getPeakArea())/(Math.log10(maxPeakArea)*iBionCurCount);//得到对应b ion位置的峰面积
//                    dMatch += Math.log10(lMSTwotrails.get(ibpos).getPeaksSum());//得到对应b ion位置的峰

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
                {
                    dMatch*=(1-minBIonDistance);//得到理论mz与实际质荷比距离

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
                {
                    dMatch += (1-minBIonDistance)* Math.log10(lMSTwotrails.get(ibposFront).getPeakArea())/(Math.log10(maxPeakArea)*iBionCurCount);//得到理论mz与实际质荷比距离，对应位置的峰面积

                }
            }
            if(ypos.size()>yMatchSize)
            {
                yipos.add(len-i-1);
                if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.Abundances)
                {
                    dMatch += Math.log10(lMSTwotrails.get(iyposFront).getPeakArea())/(Math.log10(maxPeakArea)*iYionCurCount);//得到对应y ion位置的峰面积
//                   dMatch += Math.log10(lMSTwotrails.get(iypos).getPeaksSum());//得到对应b ion位置的峰

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance)
                {
                    dMatch*=(1-minYIonDistance);//得到理论mz与实际质荷比距离

                }else if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.CombineAbundanceAndMassDistance)
                {
                    dMatch += (1-minYIonDistance)*Math.log10(lMSTwotrails.get(iyposFront).getPeakArea())/(Math.log10(maxPeakArea)*iYionCurCount);//得到理论mz与实际质荷比距离，对应位置的峰面积

                }
            }
        }
        MatchPeptdide matchPep = null;
        int bySum = bpos.size() + ypos.size();
        if(bySum>0)
        {
            if (Utils.ScoreSumMode == Utils.ScoreSumModeEnum.MassDistance) {
                dMatch = dMatch * bySum / (lMSTwotrails.size());
            }
            if(dMatch>0.0) {
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

        }
        else
            dMatch = 0.0;//未匹配到

        //       return dMatch;
        return  matchPep ;

    }
}
