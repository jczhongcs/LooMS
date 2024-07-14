package edu.uw.waterlooms.peptideMatch;

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

public class MSTwoTrailSet  implements ReadTrailFromFile, Serializable {
//    public ArrayList<MSTwoTrail> arrMSTwoTrail = new ArrayList<MSTwoTrail>();
//    public ArrayList<MSTwoTrail> arrMSTwoTrailFromRemoveMS1Trail = new ArrayList<MSTwoTrail>();

    public List<MSTwoTrail> arrMSTwoTrail = Collections.synchronizedList(new ArrayList<MSTwoTrail>());
    public List<MSTwoTrail> arrMSTwoTrailFromRemoveMS1Trail = Collections.synchronizedList(new ArrayList<MSTwoTrail>());


    public double mzHigh;
    public double mzLow;
    public double[] arrMSTwoRts;


    public Map<Double, List<MSTwoTrail>> mapRTMStwoTrail;

    public Map<Double, List<MSTwoTrail>> mapRTMStwoTrailFromRemoveMS1Trail;

    public Map<Double, Double> mapRTMSTwoMaxPeakArea;
    public int iWindowSize;

    public int iOriCount;
    public int iFinalCount;





    public void readFromeTrailFile(String trailFileName,
                                   TreeMap<Double, List<MSOneTrail>> mapRTMSoneFeature) throws IOException {

        FileReader freader;
        try {
            freader = new FileReader(trailFileName);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
            String line;

            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) { continue; }
                // if see end, add cur to list, next
                if(lineType == 0)//skip title line
                {
                    lineType = 1;
                    continue;
                }

                String[] range = line.split("\\s+");
                if(Utils.bzRemoveTrailFromSamePrecursor) {
                    //Filter out the trails in range between mzlow and mzhigh
                    //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.



                    if (Double.parseDouble(range[0]) < mzLow || Double.parseDouble(range[0]) > mzHigh) {

//                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                        MSTwoTrail cur = new MSTwoTrail(Double.parseDouble(range[0]),
                                Double.parseDouble(range[1]),
                                range[2],
                                range[3],
                                range[4],
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6])
                        );
//                        bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\t" + "peakSum" + "\t" + "peakArea" + "\n");
/*                        MSTwoTrail cur = new MSTwoTrail(
                                Double.parseDouble(range[1]),
                                Double.parseDouble(range[0]),
                                Double.parseDouble(range[2]),
                                range[4],
                                range[3],
                                range[5],
                                Double.parseDouble(range[6]),
                                Double.parseDouble(range[7])
                        );
*/
                        arrMSTwoTrail.add(cur);
                    }else
                    {
                        boolean bzMS1Trail = false;
                        if(mapRTMSoneFeature.floorKey(Double.parseDouble(range[1]))!=null) {
                            List<MSOneTrail> listRTMS1Features = mapRTMSoneFeature.get(mapRTMSoneFeature.floorKey(Double.parseDouble(range[1])));
                            int iPosBegin = Utils.binarySearch0(listRTMS1Features, 0, listRTMS1Features.size(), mzLow);//pepsMass[0]);//- Utils.H2OMass);

                            iPosBegin = Utils.getiPosBeginOfNearest(listRTMS1Features.size(), iPosBegin);//从大于等于threshold开始
                            //从小于它的最后一个开始，并考虑所有的相等的情况
                            for (int iPos = iPosBegin;
                                 iPos < listRTMS1Features.size() && listRTMS1Features.get(iPos).getMz() < mzHigh;//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                 iPos++) {
                                for (int iMassPos = 0; iMassPos < listRTMS1Features.get(iPos).getArrMS1Mzs().length; iMassPos++) {
                                    if (listRTMS1Features.get(iPos).getArrMS1Mzs()[iMassPos] == Double.parseDouble(range[0])) {
                                        bzMS1Trail = true;
                                        break;
                                    }
                                }
                                if (bzMS1Trail) break;
                            }
                        }/*else
                        {
                            System.out.println(Double.parseDouble(range[1])+"--------"+mapRTMSoneFeature.get(mapRTMSoneFeature.lastKey()).get(0).getRt());
                        }*/
                        if(!bzMS1Trail)
                        {
                            //                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                        MSTwoTrail cur = new MSTwoTrail(Double.parseDouble(range[0]),
                                Double.parseDouble(range[1]),
                                range[2],
                                range[3],
                                range[4],
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6])
                        );
//                        bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\t" + "peakSum" + "\t" + "peakArea" + "\n");
                          /*  MSTwoTrail cur = new MSTwoTrail(
                                    Double.parseDouble(range[1]),
                                    Double.parseDouble(range[0]),
                                    Double.parseDouble(range[2]),
                                    range[4],
                                    range[3],
                                    range[5],
                                    Double.parseDouble(range[6]),
                                    Double.parseDouble(range[7])
                            );

                           */
                            arrMSTwoTrail.add(cur);
                        }else//添加到移除的MS1Trail里面
                        {
                            if(Utils.bzCheckMS2TrailContainMS1Trail) {
                                //                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                        MSTwoTrail cur = new MSTwoTrail(Double.parseDouble(range[0]),
                                Double.parseDouble(range[1]),
                                range[2],
                                range[3],
                                range[4],
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6])
                        );
//                        bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\t" + "peakSum" + "\t" + "peakArea" + "\n");
                              /*  MSTwoTrail cur = new MSTwoTrail(
                                        Double.parseDouble(range[1]),
                                        Double.parseDouble(range[0]),
                                        Double.parseDouble(range[2]),
                                        range[4],
                                        range[3],
                                        range[5],
                                        Double.parseDouble(range[6]),
                                        Double.parseDouble(range[7])
                                );

                               */
                                arrMSTwoTrailFromRemoveMS1Trail.add(cur);
                            }

                        }
                    }

                }else
                {
                    //                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                        MSTwoTrail cur = new MSTwoTrail(Double.parseDouble(range[0]),
                                Double.parseDouble(range[1]),
                                range[2],
                                range[3],
                                range[4],
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6])
                        );
//                        bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\t" + "peakSum" + "\t" + "peakArea" + "\n");
                    /*MSTwoTrail cur = new MSTwoTrail(
                            Double.parseDouble(range[1]),
                            Double.parseDouble(range[0]),
                            Double.parseDouble(range[2]),
                            range[4],
                            range[3],
                            range[5],
                            Double.parseDouble(range[6]),
                            Double.parseDouble(range[7])
                    );

                     */
                    arrMSTwoTrail.add(cur);
                }

            }
            br.close();
            freader.close();
        }






//        mapRTMStwoTrail =   arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList())) ;

        //0915 用于保留只有单峰的，且具有比较高丰度的MS2 TRAIL(比当前时间trail最高面积一半要大)
        mapRTMSTwoMaxPeakArea = arrMSTwoTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.reducing(
                        0.0,
                        MSTwoTrail::getPeakArea,
                        Double::max)));


        iOriCount =arrMSTwoTrail.size();

        if(Utils.dThresholdMS2TrailPeakareaForOnePeak>0) {
            arrMSTwoTrail=arrMSTwoTrail.stream().filter(x->
                    (x.getMzs().length>1) || (x.getMzs().length==1 && x.getPeakArea()>(mapRTMSTwoMaxPeakArea.get(x.getRtApex())/Utils.dThresholdMS2TrailPeakareaForOnePeak))
            ).collect(Collectors.toCollection(ArrayList::new));
        }

        iFinalCount = arrMSTwoTrail.size();

        arrMSTwoRts = arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
                .distinct()
                .sorted()
                .toArray();

        mapRTMStwoTrail =  Collections.synchronizedMap( arrMSTwoTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList())) );

        if(Utils.bzCheckMS2TrailContainMS1Trail) {

            mapRTMStwoTrailFromRemoveMS1Trail = arrMSTwoTrailFromRemoveMS1Trail.stream().sorted()
                    .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
        }
//        mapRTMSTwoMaxPeakArea = arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.maxBy(Comparator.comparing(MSTwoTrail::getPeakArea))));







        //localHalfRank of a trail: within +/-50Da, +/- 5seconds, # of trails with peak area at least half of the current peak area.
        for(MSTwoTrail ms2trail:arrMSTwoTrail)
        {
            int ipeakAreaHalfRank = 0;
            //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
            int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, ms2trail.getRtApex() - Utils.thresholdRT);
            iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
            if (iTwoRTPosBegin > 0)
                iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position

            for (int iRTPos = iTwoRTPosBegin;
                 iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < ms2trail.getRtApex() + Utils.thresholdRT;
                 iRTPos++) {
                ipeakAreaHalfRank += mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().filter((x)->
                    x.getMzApex()>=ms2trail.getMzApex()-50.0 &&  x.getMzApex()<=ms2trail.getMzApex()+50.0 && x.getPeakArea()>=ms2trail.getPeakArea()/2
                ).count();

            }
            ms2trail.setPeakAreaHalfRank(ipeakAreaHalfRank);
        }
//        System.out.println("Test");
    }

    public void readFromeTrail(double dmz,double drt,double[] dmzs,double[] drts,double[] dints,double dpeaksSum,double dPeakArea,
                                   TreeMap<Double, List<MSOneTrail>> mapRTMSoneFeature){


        if (Utils.bzRemoveTrailFromSamePrecursor) {
            //Filter out the trails in range between mzlow and mzhigh
            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.


            if (dmz < mzLow || dmz > mzHigh) {

                MSTwoTrail cur = new MSTwoTrail(dmz,
                        drt,
                        dmzs,
                        drts,
                        dints,
                        dpeaksSum,
                        dPeakArea
                );
                arrMSTwoTrail.add(cur);
            } else {
                boolean bzMS1Trail = false;
                if (mapRTMSoneFeature.floorKey(drt) != null) {
                    List<MSOneTrail> listRTMS1Features = mapRTMSoneFeature.get(mapRTMSoneFeature.floorKey(drt));
                    int iPosBegin = Utils.binarySearch0(listRTMS1Features, 0, listRTMS1Features.size(), mzLow);//pepsMass[0]);//- Utils.H2OMass);

                    iPosBegin = Utils.getiPosBeginOfNearest(listRTMS1Features.size(), iPosBegin);//从大于等于threshold开始
                    //从小于它的最后一个开始，并考虑所有的相等的情况
                    for (int iPos = iPosBegin;
                         iPos < listRTMS1Features.size() && listRTMS1Features.get(iPos).getMz() < mzHigh;//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                         iPos++) {
                        for (int iMassPos = 0; iMassPos < listRTMS1Features.get(iPos).getArrMS1Mzs().length; iMassPos++) {
                            if (listRTMS1Features.get(iPos).getArrMS1Mzs()[iMassPos] == dmz) {
                                bzMS1Trail = true;
                                break;
                            }
                        }
                        if (bzMS1Trail) break;
                    }
                }
                if (!bzMS1Trail) {
                    //                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                    MSTwoTrail cur = new MSTwoTrail(dmz,
                            drt,
                            dmzs,
                            drts,
                            dints,
                            dpeaksSum,
                            dPeakArea
                    );
                    arrMSTwoTrail.add(cur);
                } else//添加到移除的MS1Trail里面
                {
                    if (Utils.bzCheckMS2TrailContainMS1Trail) {
                        //                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                        MSTwoTrail cur = new MSTwoTrail(dmz,
                                drt,
                                dmzs,
                                drts,
                                dints,
                                dpeaksSum,
                                dPeakArea
                        );
                        arrMSTwoTrailFromRemoveMS1Trail.add(cur);
                    }

                }
            }

        } else {

            MSTwoTrail cur = new MSTwoTrail(dmz,
                    drt,
                    dmzs,
                    drts,
                    dints,
                    dpeaksSum,
                    dPeakArea
            );
            arrMSTwoTrail.add(cur);
        }



    }

    public void ms2Process()
    {




//        mapRTMStwoTrail =   arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList())) ;

        //0915 用于保留只有单峰的，且具有比较高丰度的MS2 TRAIL(比当前时间trail最高面积一半要大)
        mapRTMSTwoMaxPeakArea = arrMSTwoTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.reducing(
                        0.0,
                        MSTwoTrail::getPeakArea,
                        Double::max)));


        iOriCount =arrMSTwoTrail.size();

        if(Utils.dThresholdMS2TrailPeakareaForOnePeak>0) {
            arrMSTwoTrail=arrMSTwoTrail.stream().filter(x->
                    (x.getMzs().length>1) || (x.getMzs().length==1 && x.getPeakArea()>(mapRTMSTwoMaxPeakArea.get(x.getRtApex())/Utils.dThresholdMS2TrailPeakareaForOnePeak))
            ).collect(Collectors.toCollection(ArrayList::new));
        }

        iFinalCount = arrMSTwoTrail.size();

        arrMSTwoRts = arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
                .distinct()
                .sorted()
                .toArray();

        mapRTMStwoTrail =  Collections.synchronizedMap( arrMSTwoTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList())) );

        if(Utils.bzCheckMS2TrailContainMS1Trail) {

            mapRTMStwoTrailFromRemoveMS1Trail = arrMSTwoTrailFromRemoveMS1Trail.stream().sorted()
                    .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
        }
//        mapRTMSTwoMaxPeakArea = arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.maxBy(Comparator.comparing(MSTwoTrail::getPeakArea))));







        //localHalfRank of a trail: within +/-50Da, +/- 5seconds, # of trails with peak area at least half of the current peak area.
        for(MSTwoTrail ms2trail:arrMSTwoTrail)
        {
            int ipeakAreaHalfRank = 0;
            //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
            int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, ms2trail.getRtApex() - Utils.thresholdRT);
            iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
            if (iTwoRTPosBegin > 0)
                iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position

            for (int iRTPos = iTwoRTPosBegin;
                 iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < ms2trail.getRtApex() + Utils.thresholdRT;
                 iRTPos++) {
                ipeakAreaHalfRank += mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().filter((x)->
                        x.getMzApex()>=ms2trail.getMzApex()-50.0 &&  x.getMzApex()<=ms2trail.getMzApex()+50.0 && x.getPeakArea()>=ms2trail.getPeakArea()/2
                ).count();

            }
            ms2trail.setPeakAreaHalfRank(ipeakAreaHalfRank);
        }
//        System.out.println("Test");
    }
    @Override
    public void readFromeTrailFile(String trailFileName) throws IOException {
        FileReader freader;
        try {
            freader = new FileReader(trailFileName);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
            String line;

            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) { continue; }
                // if see end, add cur to list, next
                if(lineType == 0)//skip title line
                {
                    lineType = 1;
                    continue;
                }

                String[] range = line.split("\\s+");
                if(Utils.bzRemoveTrailFromSamePrecursor) {
                    //Filter out the trails in range between mzlow and mzhigh
                    if (Double.parseDouble(range[0]) < mzLow || Double.parseDouble(range[0]) > mzHigh) {
                        //                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                        MSTwoTrail cur = new MSTwoTrail(Double.parseDouble(range[0]),
                                Double.parseDouble(range[1]),
                                range[2],
                                range[3],
                                range[4],
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6])
                        );
//                        bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\t" + "peakSum" + "\t" + "peakArea" + "\n");
                       /* MSTwoTrail cur = new MSTwoTrail(
                                Double.parseDouble(range[1]),
                                Double.parseDouble(range[0]),
                                Double.parseDouble(range[2]),
                                range[4],
                                range[3],
                                range[5],
                                Double.parseDouble(range[6]),
                                Double.parseDouble(range[7])
                        );

                        */
                        arrMSTwoTrail.add(cur);
                    }
                }else
                {
                    //                        String header = "mzApex\trtApex\tmzs\trts\tints\tpeaksSum\tpeakArea\tprecursorMZ\tprecursorRT\n";

                        MSTwoTrail cur = new MSTwoTrail(Double.parseDouble(range[0]),
                                Double.parseDouble(range[1]),
                                range[2],
                                range[3],
                                range[4],
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6])
                        );
//                        bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\t" + "peakSum" + "\t" + "peakArea" + "\n");
                   /* MSTwoTrail cur = new MSTwoTrail(
                            Double.parseDouble(range[1]),
                            Double.parseDouble(range[0]),
                            Double.parseDouble(range[2]),
                            range[4],
                            range[3],
                            range[5],
                            Double.parseDouble(range[6]),
                            Double.parseDouble(range[7])
                    );*/
                    arrMSTwoTrail.add(cur);
                }

            }
            br.close();
            freader.close();
        }




        //0915 用于保留只有单峰的，且具有比较高丰度的MS2 TRAIL
        mapRTMSTwoMaxPeakArea = arrMSTwoTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.reducing(
                        0.0,
                        MSTwoTrail::getPeakArea,
                        Double::max)));


        if(Utils.dThresholdMS2TrailPeakareaForOnePeak>0) {
            arrMSTwoTrail = arrMSTwoTrail.stream().filter(x ->
                    (x.getMzs().length > 1) || (x.getMzs().length == 1 && x.getPeakArea() > (mapRTMSTwoMaxPeakArea.get(x.getRtApex()) / Utils.dThresholdMS2TrailPeakareaForOnePeak))
            ).collect(Collectors.toCollection(ArrayList::new));
        }
        arrMSTwoRts = arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
                .distinct()
                .sorted()
                .toArray();



//        mapRTMSTwoMaxPeakArea = arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.reducing(
//                        0.0,
//                        MSTwoTrail::getPeakArea,
//                        Double::max)));
//
//
//        arrMSTwoRts = arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
//                .distinct()
//                .sorted()
//                .toArray();

        mapRTMStwoTrail =  arrMSTwoTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));

//        mapRTMSTwoMaxPeakArea = arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.maxBy(Comparator.comparing(MSTwoTrail::getPeakArea))));







        //localHalfRank of a trail: within +/-50Da, +/- 5seconds, # of trails with peak area at least half of the current peak area.
        for(MSTwoTrail ms2trail:arrMSTwoTrail)
        {
            int ipeakAreaHalfRank = 0;
            //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
            int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, ms2trail.getRtApex() - Utils.thresholdRT);
            iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
            if (iTwoRTPosBegin > 0)
                iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position

            for (int iRTPos = iTwoRTPosBegin;
                 iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < ms2trail.getRtApex() + Utils.thresholdRT;
                 iRTPos++) {
                ipeakAreaHalfRank += mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().filter((x)->
                        x.getMzApex()>=ms2trail.getMzApex()-50.0 &&  x.getMzApex()<=ms2trail.getMzApex()+50.0 && x.getPeakArea()>=ms2trail.getPeakArea()/2
                ).count();

            }
            ms2trail.setPeakAreaHalfRank(ipeakAreaHalfRank);
        }
//        System.out.println("Test");
    }

    @Override
    public String toString() {
        return /*"MSTwoTrailSet{" +*/
                "arrMSTwoTrail=" + arrMSTwoTrail +
                ", mzHigh=" + mzHigh +
                ", mzLow=" + mzLow +
                ", arrMSTwoRts=" + Arrays.toString(arrMSTwoRts) +
                ", mapRTMStwoTrail=" + mapRTMStwoTrail +
                ", mapRTMSTwoMaxPeakArea=" + mapRTMSTwoMaxPeakArea ;//+
//                '}';
    }

    public MSTwoTrailSet(ArrayList<MSTwoTrail> arrMSTwoTrail, double mzHigh, double mzLow, double[] arrMSTwoRts, Map<Double, List<MSTwoTrail>> mapRTMStwoTrail, Map<Double, Double> mapRTMSTwoMaxPeakArea) {
        this.arrMSTwoTrail = arrMSTwoTrail;
        this.mzHigh = mzHigh;
        this.mzLow = mzLow;
        this.arrMSTwoRts = arrMSTwoRts;
        this.mapRTMStwoTrail = mapRTMStwoTrail;
        this.mapRTMSTwoMaxPeakArea = mapRTMSTwoMaxPeakArea;
    }

    public MSTwoTrailSet()
    {
        ;
    }


}
