package edu.uw.waterlooms.peptideMatch;



import edu.uw.waterlooms.service.ParameterService;
import me.tongfei.progressbar.ProgressBar;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

public class MSOneTrailSet  implements ReadTrailFromFile{
    public ArrayList<MSOneTrail> arrMSOneTrail = new ArrayList<MSOneTrail>();
    public ArrayList<MSOneTrail> arrRawMSOneTrail = new ArrayList<MSOneTrail>();

    double[] arrRts;
    TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail;


    int iSkipScanMS1;

    double[] arrRawRts;
    TreeMap<Double, List<MSOneTrail>> mapRTMSoneRawTrail;
    public void generateRTSandMapRT()
    {
        arrRts = arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
            .distinct()
            .sorted()
            .toArray();

        mapRTMSoneTrail = arrMSOneTrail.stream().sorted()
            .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));

    }

    public void generateRawRTSandMapRT()
    {
        arrRawRts = arrRawMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                .distinct()
                .sorted()
                .toArray();

        mapRTMSoneRawTrail = arrRawMSOneTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));

    }

    public void readFromeTrailFile(String trailFileName) throws IOException {

        TreeSet<Long> treeSetMS1ID = new TreeSet<>();
        if(Utils.bzInputSameTargetMS1Precursor)//为过滤不同的ms1 PRECURSOR
        {
            FileReader freaderMS1ID;
            try {
                 freaderMS1ID = new FileReader(Config.fastaFolder+"ms1no.csv");
                try (BufferedReader br = new BufferedReader(freaderMS1ID)) {
                    String line;
                    int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

                    while ((line = br.readLine()) != null) {
                        if(lineType == 0)//skip title line
                        {
                            lineType = 1;
                            continue;
                        }
                        treeSetMS1ID.add(Long.parseLong(line));
                    }
                }

            } catch (IOException noFile) {
                throw new FileNotFoundException();
            }
        }

        FileReader freader;
        try {
            freader = new FileReader(trailFileName);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
 //           IsotopeTrail curTrail = new IsotopeTrail();
            String line;

            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) { continue; }

                if(lineType == 0)//skip title line
                {
                    lineType = 1;
                    continue;
                }


                String[] range = line.split("\\s+");
/*                MSOneTrail cur = new MSOneTrail(Long.parseLong(range[0]),
                        Double.parseDouble(range[1]),
                        Double.parseDouble(range[2]),
                        Integer.parseInt(range[3]),
                        Integer.parseInt(range[4]),
                        Double.parseDouble(range[5]),
                        Double.parseDouble(range[6]),
                        Double.parseDouble(range[7]),
                        Double.parseDouble(range[8]),
                        Double.parseDouble(range[9]),
                        Double.parseDouble(range[10]),
                        Long.parseLong(range[11]),
                        Double.parseDouble(range[12]),
                        Double.parseDouble(range[13]),
                        Double.parseDouble(range[14]));*/
                if(Utils.bzInputSameTargetMS1Precursor)//为过滤不同的ms1 PRECURSOR
                {
                    if(treeSetMS1ID.contains(Long.parseLong(range[0])))
                    {
                        MSOneTrail cur = new MSOneTrail(Long.parseLong(range[0]),
                                Double.parseDouble(range[1]),
                                Double.parseDouble(range[2]),
                                Integer.parseInt(range[3]),
                                Integer.parseInt(range[4]),
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6]),
                                Double.parseDouble(range[7]),
                                Long.parseLong(range[8]),

                                Double.parseDouble(range[9]),
                                Double.parseDouble(range[10]),
                                Double.parseDouble(range[11]),

                                Double.parseDouble(range[12]),

                                range[13],
                                range[14],
                                range[15] );

                        //select ms1 features that quality score greater than threshold
//                if(cur.getQuality_score() >= Utils.thresholdMS1QualityScore) {
//                    arrMSOneTrail.add(cur);
//                }
                        arrMSOneTrail.add(cur);//TEST NOT FILTER
                    }
                }else{
                    MSOneTrail cur = null;
                    if (Utils.bzInputConciseMSonePrecursor)
                    {
                        double dmz =  Double.parseDouble(range[1]);
                        double drt =  Double.parseDouble(range[2]);


                            cur = new MSOneTrail(
                                    Long.parseLong(range[0]),
                                    dmz,
                                    drt,
                                    Integer.parseInt(range[3]),
                                    Integer.parseInt(range[4]),
                                    Long.parseLong(range[5]),
                                    Double.parseDouble(range[6])
                            );

                    }else {
                        cur = new MSOneTrail(Long.parseLong(range[0]),
                                Double.parseDouble(range[1]),
                                Double.parseDouble(range[2]),
                                Integer.parseInt(range[3]),
                                Integer.parseInt(range[4]),
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6]),
                                Double.parseDouble(range[7]),
                                Long.parseLong(range[8]),

                                Double.parseDouble(range[9]),
                                Double.parseDouble(range[10]),
                                Double.parseDouble(range[11]),

                                Double.parseDouble(range[12]),

                                range[13],
                                range[14],
                                range[15] );
                    }


                    //select ms1 features that quality score greater than threshold
//                if(cur.getQuality_score() >= Utils.thresholdMS1QualityScore) {
//                    arrMSOneTrail.add(cur);
//                }
                    arrMSOneTrail.add(cur);//TEST NOT FILTER
                }


            }
            br.close();
            freader.close();
        }
    }
    public void readFromeTrailFile(String trailFileName,boolean bzRawFile) throws IOException {

        TreeSet<Long> treeSetMS1ID = new TreeSet<>();
        if(Utils.bzInputSameTargetMS1Precursor)//为过滤不同的ms1 PRECURSOR
        {
            FileReader freaderMS1ID;
            try {
                freaderMS1ID = new FileReader(Config.fastaFolder+"ms1no.csv");
                try (BufferedReader br = new BufferedReader(freaderMS1ID)) {
                    String line;
                    int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

                    while ((line = br.readLine()) != null) {
                        if(lineType == 0)//skip title line
                        {
                            lineType = 1;
                            continue;
                        }
                        treeSetMS1ID.add(Long.parseLong(line));
                    }
                }

            } catch (IOException noFile) {
                throw new FileNotFoundException();
            }
        }

        FileReader freader;
        try {
            freader = new FileReader(trailFileName);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
            //           IsotopeTrail curTrail = new IsotopeTrail();
            String line;

            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

            long ms1Count= arrMSOneTrail.size();
            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) { continue; }

                if(lineType == 0)//skip title line
                {
                    lineType = 1;
                    continue;
                }


                String[] range = line.split("\\s+");

                if(Utils.bzInputSameTargetMS1Precursor)//为过滤不同的ms1 PRECURSOR
                {
                    if(treeSetMS1ID.contains(Long.parseLong(range[0])))
                    {
                        MSOneTrail cur = new MSOneTrail(Long.parseLong(range[0]),
                                Double.parseDouble(range[1]),
                                Double.parseDouble(range[2]),
                                Integer.parseInt(range[3]),
                                Integer.parseInt(range[4]),
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6]),
                                Double.parseDouble(range[7]),
                                Long.parseLong(range[8]),

                                Double.parseDouble(range[9]),
                                Double.parseDouble(range[10]),
                                Double.parseDouble(range[11]),

                                Double.parseDouble(range[12]),

                                range[13],
                                range[14],
                                range[15] );

                        //select ms1 features that quality score greater than threshold
//                if(cur.getQuality_score() >= Utils.thresholdMS1QualityScore) {
//                    arrMSOneTrail.add(cur);
//                }
                        arrMSOneTrail.add(cur);//TEST NOT FILTER
                    }
                }else{
                    MSOneTrail cur = null;

                    if (bzRawFile)//(Utils.bzInputConciseMSonePrecursor)
                    {
                        double dmz =  Double.parseDouble(range[1]);
                        double drt =  Double.parseDouble(range[2]);
                        int iz = Integer.parseInt(range[3]);

                        //840.5339965820312  2773.427  53.23167
//                        if(dmz>925.482
//                                && dmz<925.483
//                        )
//                            System.out.println("debug");
                        //是否是只有原始数据
                        if(Utils.bzOnlyRawdata)
                        {
                            cur = new MSOneTrail(
                                    Long.parseLong(range[0]),
                                    dmz,
                                    drt,
                                    iz,
                                    Integer.parseInt(range[4]),
                                    Long.parseLong(range[5]),
                                    Double.parseDouble(range[6])
                            );
                            arrRawMSOneTrail.add(cur);//TEST NOT FILTER
                        }else
                        {
                            if(!isExist(dmz,drt,iz)) {

                                cur = new MSOneTrail(
                                        Long.parseLong(range[0])+ms1Count,//使precursor featureid 唯一0810
                                        dmz,
                                        drt,
                                        iz,
                                        Integer.parseInt(range[4]),
                                        Long.parseLong(range[5]),
                                        Double.parseDouble(range[6])
                                );
                                arrRawMSOneTrail.add(cur);//TEST NOT FILTER

                            }
                        }
                    }else {
                        cur = new MSOneTrail(Long.parseLong(range[0]),
                                Double.parseDouble(range[1]),
                                Double.parseDouble(range[2]),
                                Integer.parseInt(range[3]),
                                Integer.parseInt(range[4]),
                                Double.parseDouble(range[5]),
                                Double.parseDouble(range[6]),
                                Double.parseDouble(range[7]),
                                Long.parseLong(range[8]),

                                Double.parseDouble(range[9]),
                                Double.parseDouble(range[10]),
                                Double.parseDouble(range[11]),

                                Double.parseDouble(range[12]),

                                range[13],
                                range[14],
                                range[15] );
                        cur.bAddPrecursor = true;//加入的是xiangyuan数据
                        arrMSOneTrail.add(cur);//TEST NOT FILTER

                    }


                    //select ms1 features that quality score greater than threshold
//                if(cur.getQuality_score() >= Utils.thresholdMS1QualityScore) {
//                    arrMSOneTrail.add(cur);
//                }
                }


            }
            br.close();
            freader.close();
        }

        //call combine raw isotope data
        //thresholdMS1ScanSpan great than 0 then combine
        if(bzRawFile && Utils.thresholdMS1ScanSpan>0)
            combineRawIsotopeWithSkipSpan();// combineRawIsotope();
    }


    public void readFromeListTrail(List<Double[]> list,boolean bzRawFile) throws IOException {


        int mz_index = ParameterService.getMzIndex();
        int rt_index = ParameterService.getRtIndex();
        int z_index = ParameterService.getzIndex();
        int isonum_index = ParameterService.getIsonumIndex();
        int int_shape_index = ParameterService.getIntShapeIndex();
        int iso_distr_index = ParameterService.getIsoDistrIndex();
        int intensity_area_percentage_index = ParameterService.getIntensityAreaPercentageIndex();
        int scan_num_index = ParameterService.getScanSumIndex();
        int quantification_peaks_sum_index = ParameterService.getPeaksSumIndex();
        int quantification_peaks_area_index = ParameterService.getPeaksAreaIndex();
        int svr_index = ParameterService.getSvrIndex();
        int quality_index = ParameterService.getQualityIndex();
        int invalidVal = ParameterService.getInvalidVal();
        Integer id = 1;
        long ms1Count= arrMSOneTrail.size();


        TreeSet<Long> treeSetMS1ID = new TreeSet<>();
        if (Utils.bzInputSameTargetMS1Precursor)//为过滤不同的ms1 PRECURSOR
        {
            FileReader freaderMS1ID;
            try {
                freaderMS1ID = new FileReader(Config.fastaFolder + "ms1no.csv");
                try (BufferedReader br = new BufferedReader(freaderMS1ID)) {
                    String line;
                    int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

                    while ((line = br.readLine()) != null) {
                        if (lineType == 0)//skip title line
                        {
                            lineType = 1;
                            continue;
                        }
                        treeSetMS1ID.add(Long.parseLong(line));
                    }
                }

            } catch (IOException noFile) {
                throw new FileNotFoundException();
            }
        }


        for (int i = 0; i < list.size(); i++) {
            Integer z = list.get(i)[z_index].intValue();
            Integer iso_num = list.get(i)[isonum_index].intValue();
            int scannum = list.get(i)[scan_num_index].intValue();
            String mzsStr = "[" + list.get(i)[quality_index + 1].toString();
            String rtsStr = "[" + list.get(i)[quality_index + 1 + scannum].toString();
            String intsStr = "[" + list.get(i)[quality_index + 1 + scannum * 2].toString();
            for (int k = 1; k < scannum; k++) {
                mzsStr += "," + list.get(i)[quality_index + 1 + k];
                rtsStr += "," + list.get(i)[quality_index + 1 + scannum + k];
                intsStr += "," + list.get(i)[quality_index + 1 + scannum * 2 + k];
            }
            mzsStr += "]";
            rtsStr += "]";
            intsStr += "]";




            if (Utils.bzInputSameTargetMS1Precursor)//为过滤不同的ms1 PRECURSOR
            {
                if (treeSetMS1ID.contains(id)) {
//                MSOneTrail(long lid, double dmz, double drt,
//                int iz, int iisotope_num, double dintensity_shape_score,
//                double disotope_distribution_score, double dintensity_area_percentage, long lscan_num,
//                double dquantification_peaks_sum, double dquantification_peaks_area, double dsvr_score,
//                double dquality_score, String smz, String srt, String sints)

                    MSOneTrail cur = new MSOneTrail(id,
                            list.get(i)[mz_index],
                            list.get(i)[rt_index],
                            z,
                            iso_num,
                            list.get(i)[int_shape_index],
                            list.get(i)[iso_distr_index],
                            list.get(i)[intensity_area_percentage_index],
                            scannum,
                            list.get(i)[quantification_peaks_sum_index],
                            list.get(i)[quantification_peaks_area_index],
                            list.get(i)[svr_index],
                            list.get(i)[quality_index],

                            mzsStr,
                            rtsStr,
                            intsStr);
                    arrMSOneTrail.add(cur);//TEST NOT FILTER
                }
            } else {
                MSOneTrail cur = null;

                if (bzRawFile)//(Utils.bzInputConciseMSonePrecursor)
                {
                    double dmz = list.get(i)[mz_index];
                    double drt = list.get(i)[rt_index];
                    int iz = z;

                    //是否是只有原始数据
                    if (Utils.bzOnlyRawdata) {
//                        MSOneTrail(long lid, double dmz, double drt,
//                        int iz, int iisotope_num, long lscan_num,
//                        double dquantification_peaks_area)
                        cur = new MSOneTrail(
                                id,
                                dmz,
                                drt,
                                iz,
                                iso_num,
                                scannum,
                                list.get(i)[quantification_peaks_area_index]
                        );
                        arrRawMSOneTrail.add(cur);//TEST NOT FILTER
                    } else {
                        if (!isExist(dmz, drt, iz)) {

                            cur = new MSOneTrail(
                                    id + ms1Count,//使precursor featureid 唯一0810
                                    dmz,
                                    drt,
                                    iz,
                                    iso_num,
                                    scannum,
                                    list.get(i)[quantification_peaks_area_index]
                            );
                            arrRawMSOneTrail.add(cur);//TEST NOT FILTER

                        }
                    }
                } else {
                    cur = new MSOneTrail(id,
                            list.get(i)[mz_index],
                            list.get(i)[rt_index],
                            z,
                            iso_num,
                            list.get(i)[int_shape_index],
                            list.get(i)[iso_distr_index],
                            list.get(i)[intensity_area_percentage_index],
                            scannum,
                            list.get(i)[quantification_peaks_sum_index],
                            list.get(i)[quantification_peaks_area_index],
                            list.get(i)[svr_index],
                            list.get(i)[quality_index],
                            mzsStr,
                            rtsStr,
                            intsStr);
                    cur.bAddPrecursor = true;//加入的是xiangyuan数据
                    arrMSOneTrail.add(cur);//TEST NOT FILTER

                }
            }


            id++;
        }




        //call combine raw isotope data
        //thresholdMS1ScanSpan great than 0 then combine
        if (bzRawFile && Utils.thresholdMS1ScanSpan > 0)
            combineRawIsotopeWithSkipSpan();// combineRawIsotope();
    }
    public void readFromeTrail(double dmz,double drt,int iz,int iso_num,int scannum,double darea,boolean bzRawFile) {
        long ms1Count = arrMSOneTrail.size()+arrRawMSOneTrail.size();
        TreeSet<Long> treeSetMS1ID = new TreeSet<>();
        MSOneTrail cur = null;
        if (bzRawFile)//(Utils.bzInputConciseMSonePrecursor)
        {
            if (!isExist(dmz, drt, iz)) {

                cur = new MSOneTrail(
                        1 + ms1Count,//使precursor featureid 唯一0810
                        dmz,
                        drt,
                        iz,
                        iso_num,
                        scannum,
                        darea
                );
                arrRawMSOneTrail.add(cur);//TEST NOT FILTER
            }
        }
    }

    public void combinRaw() {
         //call combine raw isotope data
         //thresholdMS1ScanSpan great than 0 then combine
         if ( Utils.thresholdMS1ScanSpan > 0)
             combineRawIsotopeWithSkipSpan();// combineRawIsotope();
     }


    private void combineRawIsotopeWithSplit()
    {
        generateRawRTSandMapRT();

        ProgressBar pb = new ProgressBar("Combine Precursor Progress", arrRawRts.length);
        pb.start();
        for(int rt=0;rt<arrRawRts.length;rt++)
        {
            pb.step();
            for(MSOneTrail msoneTrail:mapRTMSoneRawTrail.get(arrRawRts[rt])) {


                if (msoneTrail.isIzCombine())
                    continue;
                double dmzMin = msoneTrail.getMz() * (1 - Utils.thresholdMS1PPM / 2.0);
                double dmzMax = msoneTrail.getMz() * (1 + Utils.thresholdMS1PPM / 2.0);
                int iz = msoneTrail.getZ();


//                if(msoneTrail.getRt()>133.4071 && msoneTrail.getRt()<133.4072
//                        && msoneTrail.getMz()>334.044 && msoneTrail.getMz()<334.045)
//                {
//
//                    int test=0;
//                }

                //initial data
                double dmaxInt = 0.0;
                ArrayList<Double> listMZs = new ArrayList<>();
                ArrayList<Double> listRTs = new ArrayList<>();
                ArrayList<Double> listINTs = new ArrayList<>();
                dmaxInt = msoneTrail.getQuantification_peaks_area();
                double dpreIntValue = msoneTrail.getQuantification_peaks_area();
                listMZs.add(msoneTrail.getMz());
                listRTs.add(msoneTrail.getRt());
                listINTs.add(msoneTrail.getQuantification_peaks_area());
                msoneTrail.setIzCombine(true);
                MSOneTrail addMSoneTrail = msoneTrail;

                int curRT =rt;
//                int findRT = rt+1;
                int findRT = rt+iSkipScanMS1;
                boolean isDecrease = false;//判断丰度是否下降，若丰度下降双上升则生成新的precursor，规则是丰度大于之前最大丰度的一倍

                ArrayList<MSOneTrail> arrRecoverCombine = new ArrayList<>();
                ArrayList<Integer> arriLowestPos = new ArrayList<>();
                int iLowestPos =0;
//                while(findRT<arrRawRts.length && findRT-curRT<= Utils.thresholdMS1ScanSpan)
                while(findRT<arrRawRts.length && findRT-curRT<= Utils.thresholdMS1ScanSpan*iSkipScanMS1)
                {
                    MSOneTrail sameMSoneTrail = canFindSameTrailInNextScan(findRT,dmzMin,dmzMax,msoneTrail.getMz(),iz);
                    if(sameMSoneTrail!=null)
                    {
                        if (sameMSoneTrail.getQuantification_peaks_area() < dpreIntValue)
                        {
                            isDecrease = true;
                            iLowestPos = findRT;

                        }else if(isDecrease)
                        {
                            if (sameMSoneTrail.getQuantification_peaks_area() > dmaxInt) {

                                if(sameMSoneTrail.getQuantification_peaks_area()*Utils.thresholdMS1IntensityPercentagofMax > dmaxInt)
                                {
                                    //split new precursor
                                    arriLowestPos.add(iLowestPos);
                                    arrRecoverCombine.add(addMSoneTrail);
                                    isDecrease = false;

                                }
                                dmaxInt = sameMSoneTrail.getQuantification_peaks_area();
                                addMSoneTrail = sameMSoneTrail;

                            }

                        }else{
                            if (sameMSoneTrail.getQuantification_peaks_area() > dmaxInt) {
                                dmaxInt = sameMSoneTrail.getQuantification_peaks_area();
                                addMSoneTrail =sameMSoneTrail;

                            }
                        }

                        listMZs.add(sameMSoneTrail.getMz());
                        listRTs.add(sameMSoneTrail.getRt());
                        listINTs.add(sameMSoneTrail.getQuantification_peaks_area());

                        sameMSoneTrail.setIzCombine(true);
                        curRT = findRT;
//                        findRT++;
                        findRT+=iSkipScanMS1;
                        dpreIntValue = sameMSoneTrail.getQuantification_peaks_area();
                    }
                    else
                    {
//                        findRT++;
                        findRT+=iSkipScanMS1;
                        //没有匹配，则为0值
                        dpreIntValue = 0;
                        isDecrease = true;
                    }

                }
                arrRecoverCombine.add(addMSoneTrail);

                if(arriLowestPos.size()<1) {
                    addMSoneTrail.setArrMS1Mzs(listMZs);
                    addMSoneTrail.setArrMS1Rts(listRTs);
                    addMSoneTrail.setArrMS1Ins(listINTs);
                    addMSoneTrail.quantifyPeaks();

                    arrMSOneTrail.add(addMSoneTrail);
                }else
                {
                    int icurr =0;
                    for(int ipos =0;ipos<=arriLowestPos.size();ipos++)
                    {
                        int ilentrail;
                        if(ipos==arriLowestPos.size())
                            ilentrail = listMZs.size();
                        else
                            ilentrail = arriLowestPos.get(ipos);

                        addMSoneTrail = arrRecoverCombine.get(ipos);

                        addMSoneTrail.setArrMS1Mzs((ArrayList<Double>) listMZs.subList(icurr,ilentrail));
                        addMSoneTrail.setArrMS1Rts((ArrayList<Double>) listRTs.subList(icurr,ilentrail));
                        addMSoneTrail.setArrMS1Ins((ArrayList<Double>) listINTs.subList(icurr,ilentrail));
                        addMSoneTrail.quantifyPeaks();
                        icurr=ilentrail;

                        arrMSOneTrail.add(addMSoneTrail);
                    }
                }

            }
        }
        pb.stop();
//        ArrayList<MSOneTrail>  sortedMSOneTrail = (ArrayList<MSOneTrail>) arrRawMSOneTrail.stream()
//                .sorted(compareByZandMZandRT)
//                .collect(Collectors.toList());
//        System.out.println("test");
    }
    /**
     * combine TRAILS based on raw signal data before ISOTOPE;
     * calculate
     */
    private void combineRawIsotope()
    {
        generateRawRTSandMapRT();
//        Comparator<MSOneTrail> compareByZandMZandRT = Comparator
//                .comparing(MSOneTrail::getRt)
//                .thenComparing(MSOneTrail::getMz);
        //filter out the ms1 precursor with mz range and z and combine, sort by time and mz
//            ArrayList<MSOneTrail>  sortedMSOneTrail = (ArrayList<MSOneTrail>) arrRawMSOneTrail.stream().filter(
//                    x->x.getMz()>=dmzMin && x.getMz()<=dmzMax && x.getZ()==iz && x.isIzCombine()==false
//            ).sorted(compareByZandMZandRT)
//                .collect(Collectors.toList());

        ProgressBar pb = new ProgressBar("Combine Precursor Progress", arrRawRts.length);
        pb.start();
        for(int rt=0;rt<arrRawRts.length;rt++)
        {
            pb.step();
            for(MSOneTrail msoneTrail:mapRTMSoneRawTrail.get(arrRawRts[rt])) {


                if (msoneTrail.isIzCombine())
                    continue;
                double dmzMin = msoneTrail.getMz() * (1 - Utils.thresholdMS1PPM / 2.0);
                double dmzMax = msoneTrail.getMz() * (1 + Utils.thresholdMS1PPM / 2.0);
                int iz = msoneTrail.getZ();


//                if(msoneTrail.getRt()>133.4071 && msoneTrail.getRt()<133.4072
//                        && msoneTrail.getMz()>334.044 && msoneTrail.getMz()<334.045)
//                {
//
//                    int test=0;
//                }

                    //initial data
                double dmaxInt = 0.0;
                ArrayList<Double> listMZs = new ArrayList<>();
                ArrayList<Double> listRTs = new ArrayList<>();
                ArrayList<Double> listINTs = new ArrayList<>();
                dmaxInt = msoneTrail.getQuantification_peaks_area();
                listMZs.add(msoneTrail.getMz());
                listRTs.add(msoneTrail.getRt());
                listINTs.add(msoneTrail.getQuantification_peaks_area());
                msoneTrail.setIzCombine(true);
                MSOneTrail addMSoneTrail = msoneTrail;

                int curRT =rt;
                int findRT = rt+1;
                while(findRT<arrRawRts.length && findRT-curRT<= Utils.thresholdMS1ScanSpan)
                {
                    MSOneTrail sameMSoneTrail = canFindSameTrailInNextScan(findRT,dmzMin,dmzMax,msoneTrail.getMz(),iz);
                    if(sameMSoneTrail!=null)
                    {
                        if (sameMSoneTrail.getQuantification_peaks_area() > dmaxInt) {
                            dmaxInt = sameMSoneTrail.getQuantification_peaks_area();
                            addMSoneTrail =sameMSoneTrail;

                        }
                        listMZs.add(sameMSoneTrail.getMz());
                        listRTs.add(sameMSoneTrail.getRt());
                        listINTs.add(sameMSoneTrail.getQuantification_peaks_area());

                        sameMSoneTrail.setIzCombine(true);
                        curRT = findRT;
                        findRT++;
                    }
                    else
                    {
                        findRT++;
                    }

                }
                addMSoneTrail.setArrMS1Mzs(listMZs);
                addMSoneTrail.setArrMS1Rts(listRTs);
                addMSoneTrail.setArrMS1Ins(listINTs);
                addMSoneTrail.quantifyPeaks();
                addMSoneTrail.setQualityScore(curRT-rt+1);//设置ms1 precorsor得分，是有多少个共同的TRAIL，以及中心位置的isotope number

                arrMSOneTrail.add(addMSoneTrail);

            }
        }
        pb.stop();
//        ArrayList<MSOneTrail>  sortedMSOneTrail = (ArrayList<MSOneTrail>) arrRawMSOneTrail.stream()
//                .sorted(compareByZandMZandRT)
//                .collect(Collectors.toList());
//        System.out.println("test");
    }



    private MSOneTrail canFindSameTrailInNextScan(int findRT, double dmzMin, double dmzMax,double dmz,int z) {
        double dMinDistance = 1000.0;
        int ireturnpos = -1;
        double drt =arrRawRts[findRT];
        int iPosBegin = Utils.binarySearch0(mapRTMSoneRawTrail.get(drt), 0, mapRTMSoneRawTrail.get(drt).size(), dmzMin);//pepsMass[0]);//- Utils.H2OMass);


        iPosBegin = Utils.getiPosBeginOfNearest(mapRTMSoneRawTrail.get(drt).size(), iPosBegin);//从大于等于threshold开始
        //从小于它的最后一个开始，并考虑所有的相等的情况
        for (int iPos = iPosBegin;
             iPos < mapRTMSoneRawTrail.get(drt).size() && mapRTMSoneRawTrail.get(drt).get(iPos).getMz() < dmzMax;//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
             iPos++) {

            if(mapRTMSoneRawTrail.get(drt).get(iPos).isIzCombine()==false && mapRTMSoneRawTrail.get(drt).get(iPos).getZ()==z
                    && Math.abs(mapRTMSoneRawTrail.get(drt).get(iPos).getMz()-dmz)<dMinDistance)
            {
                dMinDistance = Math.abs(mapRTMSoneRawTrail.get(drt).get(iPos).getMz()-dmz);
                ireturnpos = iPos;
            }
            if(mapRTMSoneRawTrail.get(drt).get(iPos).getZ()==z && Math.abs(mapRTMSoneRawTrail.get(drt).get(iPos).getMz()-dmz)>dMinDistance)
            {
                break;
            }
        }
        if(ireturnpos>-1)
        {
            return mapRTMSoneRawTrail.get(drt).get(ireturnpos);
        }
        else
        {
            return null;
        }

    }

    //判断该mz是否存在于以前的feature中
    private boolean isExist(double dmz, double drt,int iz) {
        boolean isE = false;
        double thresholdMS1RT=0.5;
        int iMS1RTPosBegin = Arrays.binarySearch(arrRts, drt - thresholdMS1RT);//Utils.thresholdMS1RT);
        iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
        if (iMS1RTPosBegin > 0)
            iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position

        for (int iRTPos = iMS1RTPosBegin;
             iRTPos < arrRts.length && arrRts[iRTPos] < drt + thresholdMS1RT;//Utils.thresholdMS1RT;
             iRTPos++) {

            List<MSOneTrail> listRTMS1Features = mapRTMSoneTrail.get(arrRts[iRTPos]);
            int iPosBegin = Utils.binarySearch0(listRTMS1Features, 0, listRTMS1Features.size(), dmz*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

            iPosBegin = Utils.getiPosBeginOfNearest(listRTMS1Features.size(), iPosBegin);//从大于等于threshold开始
            //从小于它的最后一个开始，并考虑所有的相等的情况
            for (int iPos = iPosBegin;
                 iPos < listRTMS1Features.size() && listRTMS1Features.get(iPos).getMz() < dmz*(1+Utils.thresholdPPM);//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                 iPos++) {

                if (listRTMS1Features.get(iPos).isRTMZinRTSMZS(drt,dmz,iz))
                {
                    isE =true;
                    break;
                }
            }
            if (isE)
            {
                break;
            }

        }


        return  isE;
    }


    //判断该mz是否存在于以前的feature中
    private boolean isExistWithSkipSpan(double dmz, double drt,int iz) {
        boolean isE = false;
        double thresholdMS1RT=0.5;
        int iMS1RTPosBegin = Arrays.binarySearch(arrRts, drt - thresholdMS1RT);//Utils.thresholdMS1RT);
        iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
        if (iMS1RTPosBegin > 0)
            iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position

        for (int iRTPos = iMS1RTPosBegin;
             iRTPos < arrRts.length && arrRts[iRTPos] < drt + thresholdMS1RT;//Utils.thresholdMS1RT;
             iRTPos++) {

            List<MSOneTrail> listRTMS1Features = mapRTMSoneTrail.get(arrRts[iRTPos]);
            int iPosBegin = Utils.binarySearch0(listRTMS1Features, 0, listRTMS1Features.size(), dmz*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

            iPosBegin = Utils.getiPosBeginOfNearest(listRTMS1Features.size(), iPosBegin);//从大于等于threshold开始
            //从小于它的最后一个开始，并考虑所有的相等的情况
            for (int iPos = iPosBegin;
                 iPos < listRTMS1Features.size() && listRTMS1Features.get(iPos).getMz() < dmz*(1+Utils.thresholdPPM);//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                 iPos++) {

                if (listRTMS1Features.get(iPos).isRTMZinRTSMZS(drt,dmz,iz))
                {
                    isE =true;
                    break;
                }
            }
            if (isE)
            {
                break;
            }

        }


        return  isE;
    }


    /**
     * combine TRAILS based on raw signal data before ISOTOPE;
     * calculate
     */
    private void combineRawIsotopeWithSkipSpan()
    {
        iSkipScanMS1 =  Utils.iSkipScanMS1;

        generateRawRTSandMapRT();
//        Comparator<MSOneTrail> compareByZandMZandRT = Comparator
//                .comparing(MSOneTrail::getRt)
//                .thenComparing(MSOneTrail::getMz);
        //filter out the ms1 precursor with mz range and z and combine, sort by time and mz
//            ArrayList<MSOneTrail>  sortedMSOneTrail = (ArrayList<MSOneTrail>) arrRawMSOneTrail.stream().filter(
//                    x->x.getMz()>=dmzMin && x.getMz()<=dmzMax && x.getZ()==iz && x.isIzCombine()==false
//            ).sorted(compareByZandMZandRT)
//                .collect(Collectors.toList());

        ProgressBar pb = new ProgressBar("Combine Precursor Progress", arrRawRts.length);
        pb.start();
        for(int rt=0;rt<arrRawRts.length;rt++)
        {
            pb.step();
            for(MSOneTrail msoneTrail:mapRTMSoneRawTrail.get(arrRawRts[rt])) {


                if (msoneTrail.isIzCombine())
                    continue;
                double dmzMin = msoneTrail.getMz() * (1 - Utils.thresholdMS1PPM / 2.0);
                double dmzMax = msoneTrail.getMz() * (1 + Utils.thresholdMS1PPM / 2.0);
                int iz = msoneTrail.getZ();


//                if(msoneTrail.getRt()>133.4071 && msoneTrail.getRt()<133.4072
//                        && msoneTrail.getMz()>334.044 && msoneTrail.getMz()<334.045)
//                {
//
//                    int test=0;
//                }

                //initial data
                double dmaxInt = 0.0;
                ArrayList<Double> listMZs = new ArrayList<>();
                ArrayList<Double> listRTs = new ArrayList<>();
                ArrayList<Double> listINTs = new ArrayList<>();
                dmaxInt = msoneTrail.getQuantification_peaks_area();
                listMZs.add(msoneTrail.getMz());
                listRTs.add(msoneTrail.getRt());
                listINTs.add(msoneTrail.getQuantification_peaks_area());
                msoneTrail.setIzCombine(true);
                MSOneTrail addMSoneTrail = msoneTrail;

                int curRT =rt;
//                int findRT = rt+1;
                int findRT = rt+iSkipScanMS1;
//                while(findRT<arrRawRts.length && findRT-curRT<= Utils.thresholdMS1ScanSpan)
                while(findRT<arrRawRts.length && findRT-curRT<= Utils.thresholdMS1ScanSpan*iSkipScanMS1)
                {
                    MSOneTrail sameMSoneTrail = canFindSameTrailInNextScan(findRT,dmzMin,dmzMax,msoneTrail.getMz(),iz);
                    if(sameMSoneTrail!=null)
                    {
                        if (sameMSoneTrail.getQuantification_peaks_area() > dmaxInt) {
                            dmaxInt = sameMSoneTrail.getQuantification_peaks_area();
                            addMSoneTrail =sameMSoneTrail;

                        }
                        listMZs.add(sameMSoneTrail.getMz());
                        listRTs.add(sameMSoneTrail.getRt());
                        listINTs.add(sameMSoneTrail.getQuantification_peaks_area());

                        sameMSoneTrail.setIzCombine(true);
                        curRT = findRT;
//                        findRT++;
                        findRT+=iSkipScanMS1;
                    }
                    else
                    {
//                        findRT++;
                        findRT+=iSkipScanMS1;
                    }

                }
                addMSoneTrail.setArrMS1Mzs(listMZs);
                addMSoneTrail.setArrMS1Rts(listRTs);
                addMSoneTrail.setArrMS1Ins(listINTs);
                addMSoneTrail.quantifyPeaks();
//                addMSoneTrail.setQualityScore(curRT-rt+1);//设置ms1 precorsor得分，是有多少个共同的TRAIL，以及中心位置的isotope number
                addMSoneTrail.setQualityScore((curRT-rt)/iSkipScanMS1+1);//设置ms1 precorsor得分，是有多少个共同的TRAIL，以及中心位置的isotope number

                arrMSOneTrail.add(addMSoneTrail);

            }
        }
        pb.stop();
//        ArrayList<MSOneTrail>  sortedMSOneTrail = (ArrayList<MSOneTrail>) arrRawMSOneTrail.stream()
//                .sorted(compareByZandMZandRT)
//                .collect(Collectors.toList());
//        System.out.println("test");
    }



    private void setiSkipScanMS1(OpenMzxml of ){
        if(of.num2scan!=null)
        {
            boolean bfirst = true;
            boolean bAllWindows = false;
            double dStratMZ=0.0;
            int iSkipScan= 0;
            int iWindowSize = 0;

            for (int i:of.num2scan.keySet()) {

                for (int scan2num : of.num2scan.get(i).getChildScans()) {
                    double dMZlo = of.num2scan2.get(scan2num).getPrecursor().getMzRange().getLo();
                    if(bfirst)
                    {
                        dStratMZ = dMZlo;
                        bfirst = false;
                    }
                    else
                    {
                        if(dStratMZ==dMZlo)//已经循环了一次
                        {
                            bAllWindows =true;
                        }
                    }
                    if(bAllWindows)
                        break;

                    iWindowSize++;
                }
                if (bAllWindows)
                    break;
                iSkipScan++;
            }
            iSkipScanMS1 =  iSkipScan;
        }
    }


    public void clearRTSandMapRT() {
        arrRts = null;
        mapRTMSoneTrail = null;
        arrRawRts = null;
        arrRawMSOneTrail = null;
    }
}
