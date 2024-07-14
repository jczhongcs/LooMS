package edu.uw.waterlooms.peptideMatch;



import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import me.tongfei.progressbar.ProgressBar;
import org.json.JSONWriter;


public class RunFilter {


    public void match_TrailWithMultiRT_FitRT_PeptideAllChargeCheckAtsameRT(String mgfInfile,String rawfile, String fastaInfile, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();

        // read fasta
        ArrayList<Genome> genomes = FastaFile.ReadFile(fastaInfile);
//        ArrayList<Genome> genomes = FastaFile.ReadFileFormPeptide(fastaInfile);


        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top
//        List<Peptide> peps = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
//        double[] pepsMass = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().mapToDouble(x -> x.mass).toArray();

        List<Peptide> pepsOri = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).sorted().collect(Collectors.toCollection(ArrayList::new));


        if (Utils.bRemovePeptideFromBothTargetAndDecoy) {

            TreeMap<String, List<Peptide>> mapPepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(Peptide::getComposition, TreeMap::new, Collectors.toList()));

            System.out.println(mapPepstrListPep.size());
            ArrayList<Peptide> arrPepsRemove = new ArrayList<>();
            for (String strPep : mapPepstrListPep.keySet()) {


                boolean bDecoyType = false;
                boolean bTargetType = false;
                if (mapPepstrListPep.get(strPep).size() > 1) {
                    for (Peptide pep : mapPepstrListPep.get(strPep)) {

                        boolean isDecoy =pep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }
                }
                if (!(bDecoyType && bTargetType)) {
                    arrPepsRemove.addAll(mapPepstrListPep.get(strPep));
                }

            }
            System.out.println("The number of arrPepsRemove is: " + arrPepsRemove.size());

            pepsOri = arrPepsRemove;

        }
        List<Peptide> peps = pepsOri.stream().distinct().sorted().collect(Collectors.toCollection(ArrayList::new));

        double[] pepsMass = peps.stream().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());
        pepsOri = null;


//        FileWriter fileWriterPeptideInfo = new FileWriter("peptideInfo.csv");
//        BufferedWriter bufferedWriterPeptideInfo = new BufferedWriter(fileWriterPeptideInfo);
//        bufferedWriterPeptideInfo.write( "peptide\tmass\tdMutationRate"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
//                + "\n");
//        for (Peptide pep:peps)
//        {
//            bufferedWriterPeptideInfo.write(pep.composition+'\t'+pep.mass+'\t'+pep.dMutationRate+'\n');
//        }
//
//        bufferedWriterPeptideInfo.flush();
//        bufferedWriterPeptideInfo.close();


        ArrayList<DbMatch> resLst = new ArrayList<>();




//        Set<Peptide>
        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tmsoneMZ\tmsoneZ\tmsonePeakAreaLocalRank" +
                "\tms2rt\tpeptideLength\tMS1PeptideMassError\tdPeptideMutationRate\tibyComplemetaryCount\tdCosWithPredictMSMS\tdMS1WithMS2TrailCoscinSimilarity\tMS2TrailCoscinSimilarity\tBIon\tYIon");
        for(int iw =0;iw<Utils.iMS2Charge;iw++) {
            bufferedWriter.write("\tdlogBionPeakAreaSum_C"+(iw+1)+"\tdlogYionPeakAreaSum_C"+(iw+1)+"\tdBionMassErrorSum_C"+(iw+1)+"\tdYionMassErrorSum_C"+(iw+1)+"\tiBionPeakHalfSum_C"+(iw+1)+"\tiYionPeakHalfSum_C"+(iw+1)+"\tbionMatchCount_C"+(iw+1)+"\tyionMatchCount" +
                    "_C"+(iw+1)+"\tdBionRetentionTimeErrorSum_C"+(iw+1)+"\tdYionRetentionTimeErrorSum_C"+(iw+1)+"\tiBConsective_C"+(iw+1)+"\tiYConsective_C"+(iw+1)+"\tdBionCosWithPredictMSMS_C"+(iw+1)+"\tdYionCosWithPredictMSMS_C"+(iw+1));
        }
        bufferedWriter.write( "\tarrMatchedBion\tarrMatchedYion\tadjustBYIon\tSoredIndex\twindowSize\tdBionCos60WithPredictMSMS_C1\tdCosWithPredictMSMSMatchedPredict\tprositRT" +
                "\tibAllMatchedWithloss\tibMatchNormalAndH2Oloss\tibMatchNormalAndNH3loss\tiyAllMatchedWithloss\tiyMatchNormalAndH2Oloss\tiyMatchNormalAndNH3loss\tisRawdata\tcount"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
                + "\n");



//        Map<Integer,Map<Integer,Double>> mapTriplePeptideRTScore = new TreeMap<>();

//        BitSet mapTriplePeptideRTScore = new BitSet(peps.size());
//        BitSet[] mapTriplePeptideRTScore = new BitSet[peps.size()];
//        Map<Pair<Integer,Integer>,Double> mapTriplePeptideRTScore = new HashMap<>();
//        for (MSTwoTrail msTwoTrail : MSTWOts.arrMSTwoTrail) {
//            double[] rts = msTwoTrail.getRts();
//
//            for(double drt:rts) {
//                if (uniqueValues.add(drt))
//                {
//                    List<MSTwoTrail> listMSTwoTrail = new ArrayList<>();
//                    listMSTwoTrail.add(msTwoTrail);
//
//                    mapRTMSTwoTrail.put(drt,listMSTwoTrail);
//                }
//                else
//                {
//                    mapRTMSTwoTrail.get(drt).add(msTwoTrail);
//                }
//            }
//        }

        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if(!Utils.bzOnlyRawdata)
            {
                spec.readFromeTrailFile(file,false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile,true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


//            icount = spec.arrMSOneTrail.size();
            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();

//            Map<Double, List<MSOneTrail>> mapRTMSoneTrail =  spec.arrMSOneTrail.stream().
//                    collect(Collectors.groupingBy(x-> x.getRt()));

            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));



            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if(RTWidnowDistance!=null && RTWidnowDistance.length>0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
                  IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);


//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2Intensity();
//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2IntensityFromProsit();
//            Map<String, MS2Intensity> MapMS2Intensity = new HashMap<>();//Main.getMapMS2IntensityFromProsit();//firstrunn 0829
            Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit();//firstrunn 0829

            //GET The autort information
//            Map<String, Double> MapPeptidePredictAutoRT = Main.getAutoRTInfo();
//            List<String> strpepsAutoRT = MapPeptidePredictAutoRT.keySet().stream().collect(Collectors.toList());
//            double[] dpepsAutoRT = MapPeptidePredictAutoRT.entrySet().stream().mapToDouble(x -> x.getValue()).toArray();
//
//            MapPeptidePredictAutoRT = null; //release



            //set ms1feature peakarea rank
            AtomicInteger ims1 = new AtomicInteger();
            for(MSOneTrail ms1Trail:spec.arrMSOneTrail)
            {
                long iPeakAreaRank = 1;//至少排名第一
                MSTwoTrailSet window = iswc.FindWindowWithMZ(ms1Trail.getMz(), 1);
                //get ms1feature RT distance between Math.abs(y-v) <= Utils.thresholdRT
                int iMS1RTPosBegin = Arrays.binarySearch(arrRts, ms1Trail.getRt() - Utils.thresholdMS1CheckPeakAreaRT);//正负15秒内
                iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
//                if (iMS1RTPosBegin > 0)
//                    iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position

                for (int iRTPos = iMS1RTPosBegin;
                     iRTPos < arrRts.length && arrRts[iRTPos] < ms1Trail.getRt() + Utils.thresholdMS1CheckPeakAreaRT;
                     iRTPos++) {
                    if(window!=null) {
                        //0806加入电荷相等
                        iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
                              x.getZ()==ms1Trail.getZ() &&  x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh &&
                                      x.getQuantification_peaks_area() > ms1Trail.getQuantification_peaks_area()
                        ).count();
                    }
                }
                ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));
//                if(iPeakAreaRank>0)
//                {
//                    ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));
//                }else
//                {
//                    ms1Trail.setdPeakAreaLocalRankInTimeSpan( Math.log(20.0));
//                }


/*                if(window!=null) {
                    iPeakAreaRank = spec.arrMSOneTrail.stream().filter((x) ->
                            x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh && x.getQuantification_peaks_area() >= ms1Trail.getQuantification_peaks_area()
                                    && x.getRt() >= ms1Trail.getRt() - Utils.thresholdMS1RT && x.getRt() <= ms1Trail.getRt() + Utils.thresholdMS1RT
                    ).count();
                }*/

//                System.out.println("ims1:"+ ims1++  +"----"+spec.arrMSOneTrail.size());

            }

//        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(
//                "isolationWindowRanges.txt"));
//        IsolationWindowCollection iswc = (IsolationWindowCollection) ois.readObject();

//        List<double[]> larrMSTwoRts = new ArrayList<>(iswTest.windows.size());
//        List<Map<Double, List<MSTwoTrail>>> lmapRTMStwoTrail = new ArrayList<>(iswTest.windows.size());

            // TODO: READ IN SPECTRUM RANKING FILE (Added by Caroline) Implement Rank Here


            // for each isolation window collection, match once
            //Step 1: get a msonetrail
            //Step 2: sort with RT
            //Step 3:for each msonetrail get the window
            //      check all peptide with same mass  (mz +- 30 ppm)
            //         for each mstwoTrail in window (mz in window && RT +- 5s )
            //             compare mstwotrail b y ions with peptide b y ions (mz +- 30 ppm)
            //         add peptide to msonetrail as candidate peptide with confident score
//        MSTwoTrailSet MSTWOts = new MSTwoTrailSet();
//        MSTWOts.readFromeTrailFile(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window0_ms2_trails.tsv");

//        Map<Double,List<MSTwoTrail> >mapRTMSTwoTrail = new HashMap<Double, List<MSTwoTrail>>();
            TreeSet<Double> uniqueValues = new TreeSet<>();


            //////////////////////////////////
            //for a window,initial rts and rts corresponding trails
//        double[] arrMSTwoRts = MSTWOts.arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
//                .distinct()
//                .sorted()
//                .toArray();
////        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().
////                collect(Collectors.groupingBy(x-> x.getRtApex()));
//
//
//        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
            //////////////////////////////////////

//            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = new HashMap<>();
//            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = Collections.synchronizedMap( new ConcurrentHashMap<>());
            List<MatchPeptideWithMultiCharge> listPeptideAdd = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());
            List<MatchPeptideWithMultiCharge> listPeptideAddTrue = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());


//            Map<Long, LinkedList<MatchPeptdide>> mapMSOneFeaturePeptide = new HashMap<>();

//            AtomicInteger iPrecursorNO= new AtomicInteger(1000000);
//            AtomicInteger i = new AtomicInteger();
            ProgressBar pb = new ProgressBar("Progress", arrRts.length);
            pb.start();
//            List<Double> dRTSynlist = Collections.synchronizedList(new ArrayList<>());
//            for(double dRT:arrRts)
//                dRTSynlist.add(dRT);

//            dRTSynlist.stream().parallel().forEach(
//            Arrays.stream(arrRts).forEach(
            Arrays.stream(arrRts).parallel().forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();
                        pb.step();
                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);

                        Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());


                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {
//                                        if(msonetrail.bAddPrecursor)
//                                        {


//                                        if (msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass >= pepsMass[0]) //filter m1 trail which less than the least peptide mass

                                        //begin of window with mz of m1 feature
////                                        MSTwoTrailSet window = iswc.FindWindow(msonetrail.getMass(), 2);
//                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                        // if mass matches, then pep match
//                                        if (window != null) {
//                                            double[] arrMSTwoRts = window.arrMSTwoRts;
//                                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
                                        if (msonetrail.getMass() - Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {
                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window

                                            LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide = new LinkedList<MatchPeptideWithMultiCharge>();

//                                            LinkedList<MatchPeptdide> listMaxiumPeptide = new LinkedList<MatchPeptdide>();


                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况
//                                            while (iPosBegin > 0) {
//                                                iPosBegin = iPosBegin - 1;//move to less than position
//                                                if(iPosBegin > 0 && pepsMass[iPosBegin]>pepsMass[iPosBegin-1]) break;//PEPTIDE的质量可能有相等的情况，这是peptide有序情况
//                                            }

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {



                                                //begin of window with mz of ms1 feature
                                                MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
                                                // if mass matches, then pep match
                                                MS1MatchPeptide(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msonetrail, listMaxiumPeptide, peps.get(iPos), window);
                                            }

//                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);

                                            if (listMaxiumPeptide.size() > 0)
                                            {
                                                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                                                {
                                                    listPeptideAddTrue.addAll(listMaxiumPeptide);

                                                }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
                                                {
                                                    listPeptideAddTrue.add(listMaxiumPeptide.get(0));
                                                }
                                                else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
                                                {
                                                    if (listMaxiumPeptide.size() > 1) {
                                                        listPeptideAddTrue.add(listMaxiumPeptide.get(1));
                                                    }
                                                }
                                            }
//                                            listPeptideAddTrue.addAll(listMaxiumPeptide);
                                            ims1.addAndGet(listMaxiumPeptide.size());
                                            if (Utils.bzCheckPeptideWithAllChargeAtSameRT)
                                            {
                                                //检查该MS1相同RT且已经匹配的相同PEPTIDE里面所有电荷的二级谱
                                                int iCurrenMS1Z = msonetrail.getZ();
//                                                Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());
                                                for(MatchPeptideWithMultiCharge matchPep:listMaxiumPeptide)
                                                {
                                                    //是否重复出现的PEPTIDE，若相同只进来一次
                                                    if(uniquePeptideString.contains(matchPep.pep.composition)) {
                                                        continue;
                                                    }
                                                    for(int iAllC=1;iAllC<=Utils.iMS1Charget;iAllC++)
                                                    {
                                                        //condition is the z is larger than current z
                                                        if(iAllC == iCurrenMS1Z) {
                                                            continue;//means this
                                                        }
                                                        double dPepMZWithC = Utils.MassToMz(matchPep.pep.mass,iAllC);
                                                        //search the mz in lMsoneTrail,如果存在（可能存在一个或多个）且比当前小就不search 二级谱，否则search，然后在标记该peptide已经被search，下次就不再search

                                                        if(iAllC<iCurrenMS1Z) {
                                                            //if it is matched, it should told the problem there is no matched at successive process(create a map??)
                                                            //get the window
                                                            boolean isE = false;
                                                            int iMS1PosBegin = Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            for (int iPos = iMS1PosBegin;
                                                                 iPos < lMsoneTrail.size() && lMsoneTrail.get(iPos).getMz() < dPepMZWithC*(1+Utils.thresholdPPM);//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                                 iPos++) {

                                                                isE =true;
                                                                break;

                                                            }
                                                            if(isE){
                                                                continue;
                                                            }

                                                        }else if (iAllC>iCurrenMS1Z){
                                                            //(the larger z means mz is less than current mz that has been checked before)
                                                            int iMS1PosBegin = Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            if(iMS1PosBegin
                                                                    < lMsoneTrail.size() && lMsoneTrail.get(iMS1PosBegin).getMz() < dPepMZWithC*(1+Utils.thresholdPPM))//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                            {

                                                                continue;
                                                            }


                                                        }
                                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(dPepMZWithC, 1);
                                                        MSOneTrail msOneTrail = new MSOneTrail(-1,dPepMZWithC,msonetrail.getRt(),iAllC,0,0,0);
                                                        LinkedList<MatchPeptideWithMultiCharge> listPeptide = new LinkedList<MatchPeptideWithMultiCharge>();
                                                        MS1MatchPeptide(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msOneTrail, listPeptide, matchPep.pep, window);
                                                        listPeptideAdd.addAll(listPeptide);


                                                        //check other charge



                                                    }
                                                    //check whether it has matched before

                                                    uniquePeptideString.add(matchPep.pep.composition);

                                                }


                                            }
                                        }
//                                        }


                                        //get all peps mass equals the msonetrail
//                                        peps.stream()
//                                                .filter(e->Math.abs(e.mass-msonetrail.getMass())<Utils.thresholdPPM*e.mass)
//                                                .forEach(
//                                                        pep->{
//                                                            //compare pep with all mstwotrail in rt window by using b y ions
//                                                            uniqueValues.stream().filter(v->Math.abs(y-v) <= Utils.thresholdRT).forEach(
//                                                                    twoRTtime->{
////                                                                         if(DbMatch.PepMatch(pep,mapRTMSTwoTrail.get(twoRTtime),Config.scoreMode)>0.0) {
////                                                                             System.out.println(icount.incrementAndGet() +"--------");
////                                                                             System.out.println(pep.composition + pep.mass);
////                                                                             System.out.println(msonetrail.getMass());
////                                                                             System.out.println(twoRTtime);
////                                                                         }
//
////                                                                        mapRTMSTwoTrail.get(twoRTtime).stream().forEach(
////                                                                                msTwoTrail -> {
////                                                                                    DbMatch.PepMatch(pep,msTwoTrail,Config.scoreMode);
////
////
////                                                                                }
////                                                                        );
//                                                                    }
//                                                            );
//
//
//                                                        }
//                                                );
//                                    }//end of the window with mz of ms1 feature
                                    }
                            );
                        }





//                        i.getAndIncrement();
//                        long time = System.currentTimeMillis() - start;
                        //                       System.out.println("#########time:" + i + "/" + arrRts.length);

                        //                       System.out.println(time);
                    }
            );
            pb.stop();




//            //根据RT，和预测的RT的拟合值，找出在RT范围内的PEPTIDE，添加precursor到现有的msoneTRAIL中。
//            //设置这个precursor只检查这个PEPTIDE(以前是检查所有peptide，导致非常慢，可以将rt慢长一些[-6.5,+6.5]，加1-4个电荷)
//            //考虑不同的电荷检查raw data里面是否存在误差容限内的isotope是否存在，可以检查一到两个isotope：1电荷是+1，2电荷是+0.5，3电荷是+0.33，4电荷是+0.25
//            //get  RT distance between ms1 interval 得到一个rt间隔的所有PEPTIDE
//
////            int iAutoRTPosBegin = Arrays.binarySearch(dpepsAutoRT, y-RTWidnowDistance[0]);
////            iAutoRTPosBegin = Utils.getiPosBeginOfNearest(dpepsAutoRT.length, iAutoRTPosBegin);
////            if (iAutoRTPosBegin > 0)
////                iAutoRTPosBegin = iAutoRTPosBegin - 1;//move to less than position
//            for (int iRTPos = 0;
//                 iRTPos < dpepsAutoRT.length ;
//                 iRTPos++) {
//
//                String strPep = strpepsAutoRT.get(iRTPos);
//                double dmass = Peptide.CalculateMassByPeptideStr(strPep);
//                for (int icharge = 1;icharge<=Utils.iChargeRange;icharge++)
//                {
//                    double dmz = Utils.MassToMz(dmass,icharge);
//
//                    //直接匹配相对应的PEPTIDE
//                    MSTwoTrailSet window = iswc.FindWindowWithMZ(dmz, 1);
////                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
//                    // if mass matches, then pep match
//                    if (window != null) {
//                        double[] arrMSTwoRts = window.arrMSTwoRts;
//                        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
//
//                        double rtbegin = dpepsAutoRT[iRTPos]-Utils.thresholdPredictRT;
//
//                        int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, rtbegin);
//                        iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
//                        if (iTwoRTPosBegin > 0)
//                            iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
//
//
//                        for (int iMS2RTPos = iTwoRTPosBegin;
//                             iMS2RTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < dpepsAutoRT[iRTPos] + Utils.thresholdPredictRT;
//                             iMS2RTPos++) {
//                            //只匹配一个PEPTIDE,但基于多个时间
//                            peps.get(iPos).GenerateIons();
//                            double[][] arrIntensity = null;
//                            MS2Intensity ms2Int = MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ());
//                            if (ms2Int!=null)
//                                arrIntensity = ms2Int.arrdIntensity;
//
//                            MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiCharge(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                    Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail
//                                    ,arrIntensity);
//                                                       /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
//                                                                msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
//                            if (matchPep != null) {
//                                matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错
////                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);
//
//                                if (ms2Int != null)
//                                    matchPep.dPrositRT = ms2Int.dPredictRT;
//                                matchPep.calculateCombineScore();
//                                matchPep.setDms2rt(arrMSTwoRts[iRTPos]);
//
//
//                                matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);
//                            }
//
//                        }
////                        MSOneTrail tmpPrecursor = new MSOneTrail(
////                                iPrecursorNO.getAndIncrement(), dmz, dpepsAutoRT[iRTPos],
////                                icharge, 0, 0,
////                                0, 0, 0,
////                                0, 0, 0,
////                                0, String.valueOf(dmz), String.valueOf(dpepsAutoRT[iRTPos]), "0"
////
////                        );
////                        tmpPrecursor.strGeneratePeptide = strPep;
////                        tmpPrecursor.bAddPrecursor = true;//set the precurosor bz
////                        lMsoneTrail.add(tmpPrecursor);
//                    }
//                }
//
//            }

//            for (Peptide pep: peps) {
//                // search for the window
//                IsolationWindow window = spec.FindWindow(pep.mass, 2);
//                // if mass matches, then pep match
//                if (window != null) {
//                    DbMatch res = DbMatch.PepMatch(pep, window, Config.scoreMode);
//                    if (res != null) {
//                        if (res.matchedRts.score > 0) {
//                            resLst.add(res);
//                        }
//                    }
//                }
//            }




//            List<MatchPeptdide> sorted = new LinkedList<>();
            System.out.println("add listPeptideAddTrue size:"+listPeptideAddTrue.size());

            System.out.println("add CheckPeptideWithAllChargeAtSameRT size:"+listPeptideAdd.size());
//            List<MatchPeptideWithMultiCharge> sorted = new LinkedList<>();
            List<MatchPeptideWithMultiCharge> sorted = listPeptideAdd;
            sorted.addAll(listPeptideAddTrue);
//            for(MatchPeptideWithMultiCharge map:sorted)
//            {
//                if (map==null)
//                {
//                    System.out.println("null");
//                }
//            }


//            List<MatchPeptdide> sorted = mapMSOneFeaturePeptide.entrySet().stream()
////                    .sorted(Comparator.comparing(e -> e.getValue().stream().map(MatchPeptdide::getdMatch).min(Comparator.reverseOrder()).orElse((double) 0)))
//                    //and also sort each group before collecting them in one list
//                    .flatMap(e -> e.getValue().stream().sorted(Comparator.comparing(MatchPeptdide::getdMatch))).collect(Collectors.toList());
//
//            for (MatchPeptdide mp:sorted)
//                {
//                    bufferedWriter.write(mp.toString()+"\n");
//
//                }
            //
            //
//            Set<Long> keys = mapMSOneFeaturePeptide.keySet();

//            System.out.println("matched MSonePrecursor Size:"+keys.size());
//            if (Utils.bzJustOneTrailCanMatchAPeptide) Utils.OutputMode = Utils.OutputModeEnum.OutputALL;//consider all ms1feature can match peptide, because those match score should be adjust

//            int iaddCount=0;
//            int iPrecursorCount = 0;
//            for (Long k : keys) {
//
//                List<MatchPeptideWithMultiCharge> lmp = mapMSOneFeaturePeptide.get(k);
////                List<MatchPeptdide> lmp = mapMSOneFeaturePeptide.get(k);
//                iPrecursorCount++;
//                if (lmp.size() > 0)
//                {
//                    if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
//                    {
//                        sorted.addAll(lmp);
//
//                        iaddCount+=lmp.size();
//
//                    }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
//                    {
//                        sorted.add(lmp.get(0));
//                    }
//                    else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
//                    {
//                        if (lmp.size() > 1) {
//                            sorted.add(lmp.get(1));
//                        }
//                    }
//                }
//
////                sorted.addAll(lmp);
//
//            }
//            System.out.println("precursor Count :" +iPrecursorCount);

            System.out.println("Original Count :" +ims1);
//            System.out.println("read add Count :" +iaddCount);

            System.out.println("after combine CheckPeptideWithAllChargeAtSameRT size:"+sorted.size());


           /* if(Utils.bzCalculateThePearsonCorresionWithPredictMSMS)
            {
                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致

                *//*FileReader freader;
                BufferedReader br;
                String dataFile = Config.spectrumFolder +Utils.strpDeepOutfile;//pdeep3output.txt";
                String line;
                try {
                    freader = new FileReader(dataFile);
                    br = new BufferedReader(freader);
                    while ((line = br.readLine()) != null && !(line.startsWith(">peptide|")));//定位到当前位置

                } catch (FileNotFoundException noFile) {
                    throw new FileNotFoundException();
                }*//*


                for(int iS = 0 ;iS<sorted.size();iS++) {
//                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS_JustAll(
                            MapMS2Intensity.get(
                                    sorted.get(iS).strPepcomposition+sorted.get(iS).msOneTrail.getZ()).arrdIntensity);
                    if(iS%10000==0)
                        System.out.println(iS + "---------iS---------" + sorted.get(iS).getdMatch() + sorted.get(iS).strDecoy);

                   *//* //open predict MSMS results file
                    double[][] arrPredictMSMS = new double[sorted.get(iS).pep.composition.length() - 1][4];//length -1 with b b+2 y y+2

                    String[] strPredictPepInfo = line.split("\\|");

                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序

                    while (!strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        if ((line = br.readLine()) != null && (line.startsWith(">peptide|"))) {
                            strPredictPepInfo = line.split("\\|");
                        }
                    }
                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序，过滤掉不匹配的数据

                    if (strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        String strValue = "";
                        while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) {
                            strValue = strValue + line;
                        }
                        String[] arrstrValue = strValue.replaceAll("\\[ ", "")
                                .replaceAll("\\[", "")
                                .replaceAll("\\]", "")
                                .replaceAll("  ", " ")
                                .split(" ");
                        for (int iP = 0; iP < strPredictPepInfo[1].length() - 1; iP++) {
                            for (int jP = 0; jP < 4; jP++) {
                                if(arrstrValue[iP * 4 + jP].isEmpty())
                                {
                                    System.out.println(sorted.get(iS).pep.composition);//test
                                }else
                                {
                                    arrPredictMSMS[iP][jP] = Double.parseDouble(arrstrValue[iP * 4 + jP]);

                                }
                            }
                        }
                    }
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(arrPredictMSMS);*//*
                }
               // br.close();
               // freader.close();
                //predict msms results end
            }
            System.out.println("---------CalculateThePearsonCorresionWithPredictMSMS Finish---------");
*/
            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long,Integer> mapIDCount = new HashMap();

            if (Utils.bzJustOneTrailCanMatchAPeptide) //need to adjust the match score
            {

                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致


                for(int iS = 0 ;iS<20000 && iS<sorted.size();iS++) {
                    int jS = iS + 1;
                    if (sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) {
                                l = mid + 1;
                            } else {
                                h = mid;
                            }
                        }

                        if (l - iS > 1 && sorted.get(iS).getdMatch() < sorted.get(iS + 1).getdMatch()) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS + "----jS--------------" + (jS < sorted.size() ? sorted.get(jS).getdMatch() : "last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        } else {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }
                    } else {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }
                }

                 /*   if(sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) l = mid + 1;
                            else h = mid;
                        }


//                        for(;jS<sorted.size();jS++)
//                        {
//                            if(sorted.get(iS).getdMatch() >= sorted.get(jS).getdMatch()) break;
//                        }
                        if(jS-iS>1) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS+"----jS--------------"+(jS<sorted.size()?sorted.get(jS).getdMatch():"last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        }else
                        {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }

                    }else
                    {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }

                }*/
/*                for (MatchPeptdide mp : sorted) {
                    mp.adjustScore();
                }*/


            }
            Collections.sort(sorted, Collections.reverseOrder());

            //output multiple file, concise file, only peptidefile
            FileWriter concisefileWriter;
            BufferedWriter concisebufferedWriter = null;
            FileWriter peptidesForpDeep3Writer = null;
            BufferedWriter peptidesForpDeep3bufferedWriter = null;
            FileWriter peptidesForpPrositWriter = null;
            BufferedWriter peptidesForpPrositbufferedWriter = null;
            if(Utils.bzOutputConciseAndPeptideFiles)
            {
                concisefileWriter = new FileWriter(Config.spectrumFolder + "concisefile_0502.csv");
                concisebufferedWriter = new BufferedWriter(concisefileWriter);
                peptidesForpDeep3Writer = new FileWriter(Config.spectrumFolder + "peptidesForpDeep3_0502.csv");
                peptidesForpDeep3bufferedWriter = new BufferedWriter(peptidesForpDeep3Writer);
                peptidesForpPrositWriter = new FileWriter(Config.spectrumFolder + "peptidesForpProsit_0502.csv");
                peptidesForpPrositbufferedWriter = new BufferedWriter(peptidesForpPrositWriter);
                concisebufferedWriter.write("msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tgetQuality_score\tmsonemz\tmsonecharge\tBIon\tYIon"+"\n");
                peptidesForpDeep3bufferedWriter.write("peptide\tmodinfo\tcharge"+"\n");
                peptidesForpPrositbufferedWriter.write("modified_sequence\tcollision_energy\tprecursor_charge\tfragmentation\n");

            }
            System.out.println("---------OutputALL Begin---------");

            for (MatchPeptideWithMultiCharge mp : sorted) {
                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
//                bufferedWriter.write(k+ "\t");
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");

//                    if(Utils.bzOutputConciseAndPeptideFiles)
//                    {
//
//                        concisebufferedWriter.write(mp.getConciseString()+"\n");
//                        peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
//                        peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");
//
//                    }
                }else
                {
//                    System.out.println(mp.strPepcomposition);
                    bufferedWriter.write(mp.toString() + "\n");


                }
                if(Utils.bzOutputConciseAndPeptideFiles)
                {

                    concisebufferedWriter.write(mp.getConciseString()+"\n");
                    peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
                    peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");

                }

                iCountForFDR++;

            }
            System.out.println("---------OutputALL Finish---------");

            if(Utils.bzOutputConciseAndPeptideFiles) {
                concisebufferedWriter.flush();
                peptidesForpDeep3bufferedWriter.flush();
                peptidesForpPrositbufferedWriter.flush();
                concisebufferedWriter.close();
                peptidesForpDeep3bufferedWriter.close();
                peptidesForpPrositbufferedWriter.close();
            }
        }
//        System.out.println("The number of matched peps is: " + );


        bufferedWriter.flush();
        bufferedWriter.close();

        // sort result
//        Collections.sort(resLst, Collections.reverseOrder());
//        System.out.println("The number of matched peps is: " + resLst.size());
//        bufferedWriter.flush();
//        bufferedWriter.close();

//        for (DbMatch match : resLst) {
//            // TODO Remove all System.out.println calls for production
//            System.out.println(match.id);
//            System.out.println(match.composition);
//            System.out.println(match.matchedRts.score);
//
//            // Output the information to a file
//            bufferedWriter.write(match.id + "\n");
//            bufferedWriter.write(match.composition + "\n");
//            bufferedWriter.write(Double.toString(match.matchedRts.score) + "\n");
//
//            if (Config.mode == Enums.RunMode.DEBUG) {
//                System.out.println(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")) + "\n");
//            } else {
//                System.out.println(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")) + "\n");
//            }
//        }
    }
    public void match_FDRLT01_CheckSharePeaksAtRT1MinRangeMT(String mgfInfile,String rawfile, String fastaInfile,String resultPSMFileName, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();



        //read file from top1 fdr<0.01 dataset
        ArrayList<PSMResult> arrPsmResults = SharePeaks.ReadPSMFromResultFile(resultPSMFileName);


        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);



        bufferedWriter.write("peptide\tmsonefeature\tproteinname\tdecoyTarget\tpepmass\tmsonemass\tmsonetime\tmsoneMZ" +
                "\tmsoneZ\tBIon\tYIon\ttargetscore\tcountDP\tfdr\tsharepeaksMAX\tsharePeaksAll" +
                "\tstrSharePeakMax\tstrSharePeakAll\tshareMs1feature\tshareMS1MZ\tshareMS1time\n");



        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if (!Utils.bzOnlyRawdata) {
                spec.readFromeTrailFile(file, false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile, true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();


            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));


            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if (RTWidnowDistance != null && RTWidnowDistance.length > 0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName, mapRTMSoneTrail);


            ProgressBar pb = new ProgressBar("Progress", arrPsmResults.size());
            pb.start();

            Map<PSMResult,MatchPeptideWithMultiCharge> mapPSMBYIon = Collections.synchronizedMap( new ConcurrentHashMap<>());//new HashMap<>();


            //the first psm
//            arrPsmResults.get(0).sharepeaksMax = 0;
//            MSTwoTrailSet windowTop = iswc.FindWindowWithMZ( arrPsmResults.get(0).msoneMZ, 1);
//            double[] arrMSTwoRts = windowTop.arrMSTwoRts;
//            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = windowTop.mapRTMStwoTrail;
//            MatchPeptideWithMultiCharge matchHighPSM = getBYIonsMatchWithPSM(arrPsmResults.get(0), windowTop, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
//            matchHighPSM.setMS2TrailSelected();
//            mapPSMBYIon.put(arrPsmResults.get(0),matchHighPSM);
//            bufferedWriter.write(arrPsmResults.get(0) +"\t0\t0\t0\n");

            arrPsmResults.stream().parallel().forEach(
                    y -> {
                        pb.step();
                        MSTwoTrailSet window = iswc.FindWindowWithMZ( y.msoneMZ, 1);
                        double[] arrMSTwoRts = window.arrMSTwoRts;
                        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
                        MatchPeptideWithMultiCharge matchPSM = getBYIonsMatchWithPSM(y, window, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
                        if(matchPSM==null) {
                            System.out.println(y.strpeptide+" "+y.msoneZ+" "+y.msonetime+ " "+y.msoneMZ);
                        }
                        mapPSMBYIon.put(y,matchPSM);
                    }
            );
            pb.stop();

            List<PSMResult> arrPsm = mapPSMBYIon.keySet().stream().sorted(Comparator.comparingDouble(PSMResult::getTargetscore).reversed()).collect(Collectors.toList());
            ProgressBar pbafter = new ProgressBar("Progress after", arrPsm.size());
            pbafter.start();

            for (int i = 0; i < arrPsm.size(); i++) {
                pbafter.step();
//                if(i%100000==0)
//                    System.out.println(i+"/"+arrPsm.size());
                //get rt range[-1min,1min] and the same mz window from [0,i-1]
                PSMResult curPsm = arrPsmResults.get(i);

                int iSharePeak = 0;
                double maxShareMS1MZ = 0.0;
                double maxShareMS1Time = 0.0;
                String strshareMax="";
                Long ims1id =0L;
                for (int is = 0; is < i; is++) {

                    PSMResult highPsm = arrPsmResults.get(is);
                    if (curPsm.msonetime > highPsm.msonetime - 1 && curPsm.msonetime < highPsm.msonetime + 1
                            && iswc.isSmaeMS1Window(curPsm.msoneMZ, highPsm.msoneMZ)) {

//                        MSTwoTrailSet window = iswc.FindWindowWithMZ(curPsm.msoneMZ, 1);

//                        int iCurSharePeak = getSharePeak(curPsm, highPsm,window,mapRTMSoneTrail,mapPSMBYIon,arrRts);
                        int iCurSharePeak = getSharePeak(curPsm, highPsm,mapPSMBYIon);
                        if (iCurSharePeak > iSharePeak)
                        {
                            iSharePeak = iCurSharePeak;
                            maxShareMS1MZ = highPsm.msoneMZ;
                            maxShareMS1Time = highPsm.msonetime;
                            ims1id = highPsm.msonefeature;
                            strshareMax = mapPSMBYIon.get(curPsm).strSharePeaksMaxInfo;
                        }
                    }


                }
                curPsm.sharepeaksMax = iSharePeak;
                curPsm.strSharePeaksMaxInfo = strshareMax;
                if(mapPSMBYIon.get(curPsm)!=null) {
//                    curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALL();
                    curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALLSpeed();
                    curPsm.strSharePeaksAllInfo = mapPSMBYIon.get(curPsm).strSharePeaksAllInfo;
                    mapPSMBYIon.get(curPsm).setMS2TrailSelectedSpeed();
//                    mapPSMBYIon.get(curPsm).setMS2TrailSelected();


                }
//                else {
//                    MSTwoTrailSet windowCurPSM = iswc.FindWindowWithMZ( curPsm.msoneMZ, 1);
//                    double[] arrMSTwoRtsCurPSM = windowCurPSM.arrMSTwoRts;
//                    Map<Double, List<MSTwoTrail>> mapRTMStwoTrailCurPSM = windowCurPSM.mapRTMStwoTrail;
//                    MatchPeptideWithMultiCharge matchHighPSMCurPSM = getBYIonsMatchWithPSM(curPsm, windowCurPSM, mapRTMSoneTrail, arrMSTwoRtsCurPSM, mapRTMStwoTrailCurPSM,arrRts);
//                    matchHighPSMCurPSM.setMS2TrailSelected();
//                    mapPSMBYIon.put(curPsm,matchHighPSMCurPSM);
//                }

                bufferedWriter.write(curPsm +"\t"+ims1id+"\t"+maxShareMS1MZ+"\t"+maxShareMS1Time+"\n");


            }
            pbafter.stop();


        }


        bufferedWriter.flush();
        bufferedWriter.close();

    }


    public void match_TrailWithMultiRT_FitRT_PeptideAllChargeCheckAtsameRT_CorssWindow_topPredictIntensity(String mgfInfile,String rawfile, String fastaInfile, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();

        // read fasta
        ArrayList<Genome> genomes = FastaFile.ReadFile(fastaInfile);
//        ArrayList<Genome> genomes = FastaFile.ReadFileFormPeptide(fastaInfile);


        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top
//        List<Peptide> peps = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
//        double[] pepsMass = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().mapToDouble(x -> x.mass).toArray();

        List<Peptide> pepsOri = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).sorted().collect(Collectors.toCollection(ArrayList::new));


        if (Utils.bRemovePeptideFromBothTargetAndDecoy) {

            TreeMap<String, List<Peptide>> mapPepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(Peptide::getComposition, TreeMap::new, Collectors.toList()));

            TreeMap<String, List<Peptide>> mapILreplacePepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(Peptide::getReplaceILComposition, TreeMap::new, Collectors.toList()));

            System.out.println(mapPepstrListPep.size());
            ArrayList<Peptide> arrPepsRemove = new ArrayList<>();

            for (String strPep : mapPepstrListPep.keySet()) {


                boolean bDecoyType = false;
                boolean bTargetType = false;

                if (mapPepstrListPep.get(strPep).size() > 1) {
                    for (Peptide pep : mapPepstrListPep.get(strPep)) {

                        boolean isDecoy =pep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }
                }
                //REMOVE I L MUTATE peptide in both decoy and target
                if(strPep.contains("I")||strPep.contains("L"))
                {
                    for(Peptide muSamePep:mapILreplacePepstrListPep.get(strPep.replaceAll("L","I")))
                    {
                        boolean isDecoy =muSamePep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }

                }


                if (!(bDecoyType && bTargetType)) {
                    arrPepsRemove.addAll(mapPepstrListPep.get(strPep));
                }

            }
            System.out.println("The number of arrPepsRemove is: " + arrPepsRemove.size());

            pepsOri = arrPepsRemove;

        }
        List<Peptide> peps = pepsOri.stream().distinct().sorted().collect(Collectors.toCollection(ArrayList::new));



        double[] pepsMass = peps.stream().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());
        pepsOri = null;


//        FileWriter fileWriterPeptideInfo = new FileWriter("peptideInfo.csv");
//        BufferedWriter bufferedWriterPeptideInfo = new BufferedWriter(fileWriterPeptideInfo);
//        bufferedWriterPeptideInfo.write( "peptide\tmass\tdMutationRate"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
//                + "\n");
//        for (Peptide pep:peps)
//        {
//            bufferedWriterPeptideInfo.write(pep.composition+'\t'+pep.mass+'\t'+pep.dMutationRate+'\n');
//        }
//
//        bufferedWriterPeptideInfo.flush();
//        bufferedWriterPeptideInfo.close();


        ArrayList<DbMatch> resLst = new ArrayList<>();




//        Set<Peptide>
        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tmsoneMZ\tmsoneZ\tmsonePeakAreaLocalRank" +
                "\tms2rt\tpeptideLength\tMS1PeptideMassError\tdPeptideMutationRate\tibyComplemetaryCount\tdCosWithPredictMSMS\tdMS1WithMS2TrailCoscinSimilarity\tMS2TrailCoscinSimilarity\tBIon\tYIon");
        for(int iw =0;iw<Utils.iMS2Charge;iw++) {
            bufferedWriter.write("\tdlogBionPeakAreaSum_C"+(iw+1)+"\tdlogYionPeakAreaSum_C"+(iw+1)+"\tdBionMassErrorSum_C"+(iw+1)+"\tdYionMassErrorSum_C"+(iw+1)+"\tiBionPeakHalfSum_C"+(iw+1)+"\tiYionPeakHalfSum_C"+(iw+1)+"\tbionMatchCount_C"+(iw+1)+"\tyionMatchCount" +
                    "_C"+(iw+1)+"\tdBionRetentionTimeErrorSum_C"+(iw+1)+"\tdYionRetentionTimeErrorSum_C"+(iw+1)+"\tiBConsective_C"+(iw+1)+"\tiYConsective_C"+(iw+1)+"\tdBionCosWithPredictMSMS_C"+(iw+1)+"\tdYionCosWithPredictMSMS_C"+(iw+1));
        }
        bufferedWriter.write( "\tarrMatchedBion\tarrMatchedYion\tadjustBYIon\tSoredIndex\twindowSize\tdBionCos60WithPredictMSMS_C1\tdCosWithPredictMSMSMatchedPredict\tprositRT" +
                "\tibAllMatchedWithloss\tibMatchNormalAndH2Oloss\tibMatchNormalAndNH3loss\tiyAllMatchedWithloss\tiyMatchNormalAndH2Oloss\tiyMatchNormalAndNH3loss\tisRawdata" +
                "\tpeptideAnomiAcidFreRate\tlastIonNumBC1\tlastIonNumYC1\ttop6allCos\ttop6c1Cos\ttop6Bc1Cos\ttop6Yc1Cos\tcount"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
                + "\n");



//        Map<Integer,Map<Integer,Double>> mapTriplePeptideRTScore = new TreeMap<>();

//        BitSet mapTriplePeptideRTScore = new BitSet(peps.size());
//        BitSet[] mapTriplePeptideRTScore = new BitSet[peps.size()];
//        Map<Pair<Integer,Integer>,Double> mapTriplePeptideRTScore = new HashMap<>();
//        for (MSTwoTrail msTwoTrail : MSTWOts.arrMSTwoTrail) {
//            double[] rts = msTwoTrail.getRts();
//
//            for(double drt:rts) {
//                if (uniqueValues.add(drt))
//                {
//                    List<MSTwoTrail> listMSTwoTrail = new ArrayList<>();
//                    listMSTwoTrail.add(msTwoTrail);
//
//                    mapRTMSTwoTrail.put(drt,listMSTwoTrail);
//                }
//                else
//                {
//                    mapRTMSTwoTrail.get(drt).add(msTwoTrail);
//                }
//            }
//        }

        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if(!Utils.bzOnlyRawdata)
            {
                spec.readFromeTrailFile(file,false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile,true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


//            icount = spec.arrMSOneTrail.size();
            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();

//            Map<Double, List<MSOneTrail>> mapRTMSoneTrail =  spec.arrMSOneTrail.stream().
//                    collect(Collectors.groupingBy(x-> x.getRt()));

            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));



            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if(RTWidnowDistance!=null && RTWidnowDistance.length>0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
//            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection();
            iswc.IsolationWindowCollection_paralle(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);

//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2Intensity();
//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2IntensityFromProsit();
//            Map<String, MS2Intensity> MapMS2Intensity = new HashMap<>();//Main.getMapMS2IntensityFromProsit();//firstrunn 0829
            Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit();//firstrunn 0829

            //GET The autort information
//            Map<String, Double> MapPeptidePredictAutoRT = Main.getAutoRTInfo();
//            List<String> strpepsAutoRT = MapPeptidePredictAutoRT.keySet().stream().collect(Collectors.toList());
//            double[] dpepsAutoRT = MapPeptidePredictAutoRT.entrySet().stream().mapToDouble(x -> x.getValue()).toArray();
//
//            MapPeptidePredictAutoRT = null; //release



            //set ms1feature peakarea rank
            AtomicInteger ims1 = new AtomicInteger();


            ProgressBar pbMS1 = new ProgressBar("Progress MS1 area rank", arrRts.length);
            pbMS1.start();
//            for(MSOneTrail ms1Trail:spec.arrMSOneTrail)
            spec.arrMSOneTrail.stream().parallel().forEach(ms1->
            {
                MSOneTrail ms1Trail = ms1;
                pbMS1.step();
                long iPeakAreaRank = 1;//至少排名第一
                List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(ms1Trail.getMz(), 1,true);
                for(MSTwoTrailSet  window:arrwindow) {
                    //get ms1feature RT distance between Math.abs(y-v) <= Utils.thresholdRT
                    int iMS1RTPosBegin = Arrays.binarySearch(arrRts, ms1Trail.getRt() - Utils.thresholdMS1CheckPeakAreaRT);//正负15秒内
                    iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
//                if (iMS1RTPosBegin > 0)
//                    iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position

                    for (int iRTPos = iMS1RTPosBegin;
                         iRTPos < arrRts.length && arrRts[iRTPos] < ms1Trail.getRt() + Utils.thresholdMS1CheckPeakAreaRT;
                         iRTPos++) {
                        if (window != null) {
                            //0806加入电荷相等
//                            iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
//                                    x.getZ() == ms1Trail.getZ() && x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh &&
//                                            x.getQuantification_peaks_area() > ms1Trail.getQuantification_peaks_area()
//                            ).count();
                            //0914 multi thread
                            iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
                                    x.getZ() == ms1Trail.getZ() && x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh &&
                                            x.getQuantification_peaks_area() > ms1Trail.getQuantification_peaks_area()
                            ).count();
                        }
                    }
                }
                ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));
//                if(iPeakAreaRank>0)
//                {
//                    ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));
//                }else
//                {
//                    ms1Trail.setdPeakAreaLocalRankInTimeSpan( Math.log(20.0));
//                }


/*                if(window!=null) {
                    iPeakAreaRank = spec.arrMSOneTrail.stream().filter((x) ->
                            x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh && x.getQuantification_peaks_area() >= ms1Trail.getQuantification_peaks_area()
                                    && x.getRt() >= ms1Trail.getRt() - Utils.thresholdMS1RT && x.getRt() <= ms1Trail.getRt() + Utils.thresholdMS1RT
                    ).count();
                }*/

//                System.out.println("ims1:"+ ims1++  +"----"+spec.arrMSOneTrail.size());



            });
            pbMS1.stop();

//        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(
//                "isolationWindowRanges.txt"));
//        IsolationWindowCollection iswc = (IsolationWindowCollection) ois.readObject();

//        List<double[]> larrMSTwoRts = new ArrayList<>(iswTest.windows.size());
//        List<Map<Double, List<MSTwoTrail>>> lmapRTMStwoTrail = new ArrayList<>(iswTest.windows.size());

            // TODO: READ IN SPECTRUM RANKING FILE (Added by Caroline) Implement Rank Here


            // for each isolation window collection, match once
            //Step 1: get a msonetrail
            //Step 2: sort with RT
            //Step 3:for each msonetrail get the window
            //      check all peptide with same mass  (mz +- 30 ppm)
            //         for each mstwoTrail in window (mz in window && RT +- 5s )
            //             compare mstwotrail b y ions with peptide b y ions (mz +- 30 ppm)
            //         add peptide to msonetrail as candidate peptide with confident score
//        MSTwoTrailSet MSTWOts = new MSTwoTrailSet();
//        MSTWOts.readFromeTrailFile(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window0_ms2_trails.tsv");

//        Map<Double,List<MSTwoTrail> >mapRTMSTwoTrail = new HashMap<Double, List<MSTwoTrail>>();
            TreeSet<Double> uniqueValues = new TreeSet<>();


            //////////////////////////////////
            //for a window,initial rts and rts corresponding trails
//        double[] arrMSTwoRts = MSTWOts.arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
//                .distinct()
//                .sorted()
//                .toArray();
////        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().
////                collect(Collectors.groupingBy(x-> x.getRtApex()));
//
//
//        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
            //////////////////////////////////////

//            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = new HashMap<>();
//            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = Collections.synchronizedMap( new ConcurrentHashMap<>());
            List<MatchPeptideWithMultiCharge> listPeptideAdd = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());
            List<MatchPeptideWithMultiCharge> listPeptideAddTrue = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());


//            Map<Long, LinkedList<MatchPeptdide>> mapMSOneFeaturePeptide = new HashMap<>();

//            AtomicInteger iPrecursorNO= new AtomicInteger(1000000);
//            AtomicInteger i = new AtomicInteger();
            ProgressBar pb = new ProgressBar("Progress", arrRts.length);
            pb.start();
//            List<Double> dRTSynlist = Collections.synchronizedList(new ArrayList<>());
//            for(double dRT:arrRts)
//                dRTSynlist.add(dRT);

//            dRTSynlist.stream().parallel().forEach(
//            Arrays.stream(arrRts).forEach(
            Arrays.stream(arrRts).parallel().forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();
                        pb.step();
                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);

                        Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());


                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {
//                                        if(msonetrail.bAddPrecursor)
//                                        {


//                                        if (msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass >= pepsMass[0]) //filter m1 trail which less than the least peptide mass

                                        //begin of window with mz of m1 feature
////                                        MSTwoTrailSet window = iswc.FindWindow(msonetrail.getMass(), 2);
//                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                        // if mass matches, then pep match
//                                        if (window != null) {
//                                            double[] arrMSTwoRts = window.arrMSTwoRts;
//                                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
                                        if (msonetrail.getMass() - Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {
                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window

                                            LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide = new LinkedList<MatchPeptideWithMultiCharge>();

//                                            LinkedList<MatchPeptdide> listMaxiumPeptide = new LinkedList<MatchPeptdide>();


                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况
//                                            while (iPosBegin > 0) {
//                                                iPosBegin = iPosBegin - 1;//move to less than position
//                                                if(iPosBegin > 0 && pepsMass[iPosBegin]>pepsMass[iPosBegin-1]) break;//PEPTIDE的质量可能有相等的情况，这是peptide有序情况
//                                            }

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {



                                                //begin of window with mz of ms1 feature
//                                                MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
                                                // if mass matches, then pep match
                                                List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(msonetrail.getMz(), 1,true);
                                                for(MSTwoTrailSet  window:arrwindow) {
                                                    MS1MatchPeptide_top6PredictIntensity(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msonetrail, listMaxiumPeptide, peps.get(iPos), window);
                                                }
                                            }

//                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);

                                            if (listMaxiumPeptide.size() > 0)
                                            {
                                                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                                                {
                                                    listPeptideAddTrue.addAll(listMaxiumPeptide);

                                                }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
                                                {
                                                    listPeptideAddTrue.add(listMaxiumPeptide.get(0));
                                                }
                                                else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
                                                {
                                                    if (listMaxiumPeptide.size() > 1) {
                                                        listPeptideAddTrue.add(listMaxiumPeptide.get(1));
                                                    }
                                                }
                                            }
//                                            listPeptideAddTrue.addAll(listMaxiumPeptide);
                                            ims1.addAndGet(listMaxiumPeptide.size());
                                            if (Utils.bzCheckPeptideWithAllChargeAtSameRT)
                                            {
                                                //检查该MS1相同RT且已经匹配的相同PEPTIDE里面所有电荷的二级谱
                                                int iCurrenMS1Z = msonetrail.getZ();
//                                                Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());
                                                for(MatchPeptideWithMultiCharge matchPep:listMaxiumPeptide)
                                                {
                                                    //是否重复出现的PEPTIDE，若相同只进来一次
                                                    if(uniquePeptideString.contains(matchPep.pep.composition)) {
                                                        continue;
                                                    }
                                                    for(int iAllC=1;iAllC<=Utils.iMS1Charget;iAllC++)
                                                    {
                                                        //condition is the z is larger than current z
                                                        if(iAllC == iCurrenMS1Z) {
                                                            continue;//means this
                                                        }
                                                        double dPepMZWithC = Utils.MassToMz(matchPep.pep.mass,iAllC);
                                                        //search the mz in lMsoneTrail,如果存在（可能存在一个或多个）且比当前小就不search 二级谱，否则search，然后在标记该peptide已经被search，下次就不再search

                                                        if(iAllC<iCurrenMS1Z) {
                                                            //if it is matched, it should told the problem there is no matched at successive process(create a map??)
                                                            //get the window
                                                            boolean isE = false;
                                                            int iMS1PosBegin = Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            for (int iPos = iMS1PosBegin;
                                                                 iPos < lMsoneTrail.size() && lMsoneTrail.get(iPos).getMz() < dPepMZWithC*(1+Utils.thresholdPPM);//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                                 iPos++) {

                                                                isE =true;
                                                                break;

                                                            }
                                                            if(isE){
                                                                continue;
                                                            }

                                                        }else if (iAllC>iCurrenMS1Z){
                                                            //(the larger z means mz is less than current mz that has been checked before)
                                                            int iMS1PosBegin = Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            if(iMS1PosBegin
                                                                    < lMsoneTrail.size() && lMsoneTrail.get(iMS1PosBegin).getMz() < dPepMZWithC*(1+Utils.thresholdPPM))//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                            {

                                                                continue;
                                                            }


                                                        }
//                                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(dPepMZWithC, 1);
                                                        List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(dPepMZWithC, 1,true);
                                                        for(MSTwoTrailSet  window:arrwindow) {
                                                            MSOneTrail msOneTrail = new MSOneTrail(-1, dPepMZWithC, msonetrail.getRt(), iAllC, 0, 0, 0);
                                                            LinkedList<MatchPeptideWithMultiCharge> listPeptide = new LinkedList<MatchPeptideWithMultiCharge>();
                                                            MS1MatchPeptide(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msOneTrail, listPeptide, matchPep.pep, window);
                                                            listPeptideAdd.addAll(listPeptide);
                                                        }

                                                        //check other charge



                                                    }
                                                    //check whether it has matched before

                                                    uniquePeptideString.add(matchPep.pep.composition);

                                                }


                                            }
                                        }
//                                        }


                                        //get all peps mass equals the msonetrail
//                                        peps.stream()
//                                                .filter(e->Math.abs(e.mass-msonetrail.getMass())<Utils.thresholdPPM*e.mass)
//                                                .forEach(
//                                                        pep->{
//                                                            //compare pep with all mstwotrail in rt window by using b y ions
//                                                            uniqueValues.stream().filter(v->Math.abs(y-v) <= Utils.thresholdRT).forEach(
//                                                                    twoRTtime->{
////                                                                         if(DbMatch.PepMatch(pep,mapRTMSTwoTrail.get(twoRTtime),Config.scoreMode)>0.0) {
////                                                                             System.out.println(icount.incrementAndGet() +"--------");
////                                                                             System.out.println(pep.composition + pep.mass);
////                                                                             System.out.println(msonetrail.getMass());
////                                                                             System.out.println(twoRTtime);
////                                                                         }
//
////                                                                        mapRTMSTwoTrail.get(twoRTtime).stream().forEach(
////                                                                                msTwoTrail -> {
////                                                                                    DbMatch.PepMatch(pep,msTwoTrail,Config.scoreMode);
////
////
////                                                                                }
////                                                                        );
//                                                                    }
//                                                            );
//
//
//                                                        }
//                                                );
//                                    }//end of the window with mz of ms1 feature
                                    }
                            );
                        }





//                        i.getAndIncrement();
//                        long time = System.currentTimeMillis() - start;
                        //                       System.out.println("#########time:" + i + "/" + arrRts.length);

                        //                       System.out.println(time);
                    }
            );
            pb.stop();




//            //根据RT，和预测的RT的拟合值，找出在RT范围内的PEPTIDE，添加precursor到现有的msoneTRAIL中。
//            //设置这个precursor只检查这个PEPTIDE(以前是检查所有peptide，导致非常慢，可以将rt慢长一些[-6.5,+6.5]，加1-4个电荷)
//            //考虑不同的电荷检查raw data里面是否存在误差容限内的isotope是否存在，可以检查一到两个isotope：1电荷是+1，2电荷是+0.5，3电荷是+0.33，4电荷是+0.25
//            //get  RT distance between ms1 interval 得到一个rt间隔的所有PEPTIDE
//
////            int iAutoRTPosBegin = Arrays.binarySearch(dpepsAutoRT, y-RTWidnowDistance[0]);
////            iAutoRTPosBegin = Utils.getiPosBeginOfNearest(dpepsAutoRT.length, iAutoRTPosBegin);
////            if (iAutoRTPosBegin > 0)
////                iAutoRTPosBegin = iAutoRTPosBegin - 1;//move to less than position
//            for (int iRTPos = 0;
//                 iRTPos < dpepsAutoRT.length ;
//                 iRTPos++) {
//
//                String strPep = strpepsAutoRT.get(iRTPos);
//                double dmass = Peptide.CalculateMassByPeptideStr(strPep);
//                for (int icharge = 1;icharge<=Utils.iChargeRange;icharge++)
//                {
//                    double dmz = Utils.MassToMz(dmass,icharge);
//
//                    //直接匹配相对应的PEPTIDE
//                    MSTwoTrailSet window = iswc.FindWindowWithMZ(dmz, 1);
////                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
//                    // if mass matches, then pep match
//                    if (window != null) {
//                        double[] arrMSTwoRts = window.arrMSTwoRts;
//                        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
//
//                        double rtbegin = dpepsAutoRT[iRTPos]-Utils.thresholdPredictRT;
//
//                        int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, rtbegin);
//                        iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
//                        if (iTwoRTPosBegin > 0)
//                            iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
//
//
//                        for (int iMS2RTPos = iTwoRTPosBegin;
//                             iMS2RTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < dpepsAutoRT[iRTPos] + Utils.thresholdPredictRT;
//                             iMS2RTPos++) {
//                            //只匹配一个PEPTIDE,但基于多个时间
//                            peps.get(iPos).GenerateIons();
//                            double[][] arrIntensity = null;
//                            MS2Intensity ms2Int = MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ());
//                            if (ms2Int!=null)
//                                arrIntensity = ms2Int.arrdIntensity;
//
//                            MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiCharge(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                    Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail
//                                    ,arrIntensity);
//                                                       /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
//                                                                msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
//                            if (matchPep != null) {
//                                matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错
////                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);
//
//                                if (ms2Int != null)
//                                    matchPep.dPrositRT = ms2Int.dPredictRT;
//                                matchPep.calculateCombineScore();
//                                matchPep.setDms2rt(arrMSTwoRts[iRTPos]);
//
//
//                                matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);
//                            }
//
//                        }
////                        MSOneTrail tmpPrecursor = new MSOneTrail(
////                                iPrecursorNO.getAndIncrement(), dmz, dpepsAutoRT[iRTPos],
////                                icharge, 0, 0,
////                                0, 0, 0,
////                                0, 0, 0,
////                                0, String.valueOf(dmz), String.valueOf(dpepsAutoRT[iRTPos]), "0"
////
////                        );
////                        tmpPrecursor.strGeneratePeptide = strPep;
////                        tmpPrecursor.bAddPrecursor = true;//set the precurosor bz
////                        lMsoneTrail.add(tmpPrecursor);
//                    }
//                }
//
//            }

//            for (Peptide pep: peps) {
//                // search for the window
//                IsolationWindow window = spec.FindWindow(pep.mass, 2);
//                // if mass matches, then pep match
//                if (window != null) {
//                    DbMatch res = DbMatch.PepMatch(pep, window, Config.scoreMode);
//                    if (res != null) {
//                        if (res.matchedRts.score > 0) {
//                            resLst.add(res);
//                        }
//                    }
//                }
//            }




//            List<MatchPeptdide> sorted = new LinkedList<>();
            System.out.println("add listPeptideAddTrue size:"+listPeptideAddTrue.size());

            System.out.println("add CheckPeptideWithAllChargeAtSameRT size:"+listPeptideAdd.size());
//            List<MatchPeptideWithMultiCharge> sorted = new LinkedList<>();
            List<MatchPeptideWithMultiCharge> sorted = listPeptideAdd;
            sorted.addAll(listPeptideAddTrue);
//            for(MatchPeptideWithMultiCharge map:sorted)
//            {
//                if (map==null)
//                {
//                    System.out.println("null");
//                }
//            }


//            List<MatchPeptdide> sorted = mapMSOneFeaturePeptide.entrySet().stream()
////                    .sorted(Comparator.comparing(e -> e.getValue().stream().map(MatchPeptdide::getdMatch).min(Comparator.reverseOrder()).orElse((double) 0)))
//                    //and also sort each group before collecting them in one list
//                    .flatMap(e -> e.getValue().stream().sorted(Comparator.comparing(MatchPeptdide::getdMatch))).collect(Collectors.toList());
//
//            for (MatchPeptdide mp:sorted)
//                {
//                    bufferedWriter.write(mp.toString()+"\n");
//
//                }
            //
            //
//            Set<Long> keys = mapMSOneFeaturePeptide.keySet();

//            System.out.println("matched MSonePrecursor Size:"+keys.size());
//            if (Utils.bzJustOneTrailCanMatchAPeptide) Utils.OutputMode = Utils.OutputModeEnum.OutputALL;//consider all ms1feature can match peptide, because those match score should be adjust

//            int iaddCount=0;
//            int iPrecursorCount = 0;
//            for (Long k : keys) {
//
//                List<MatchPeptideWithMultiCharge> lmp = mapMSOneFeaturePeptide.get(k);
////                List<MatchPeptdide> lmp = mapMSOneFeaturePeptide.get(k);
//                iPrecursorCount++;
//                if (lmp.size() > 0)
//                {
//                    if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
//                    {
//                        sorted.addAll(lmp);
//
//                        iaddCount+=lmp.size();
//
//                    }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
//                    {
//                        sorted.add(lmp.get(0));
//                    }
//                    else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
//                    {
//                        if (lmp.size() > 1) {
//                            sorted.add(lmp.get(1));
//                        }
//                    }
//                }
//
////                sorted.addAll(lmp);
//
//            }
//            System.out.println("precursor Count :" +iPrecursorCount);

            System.out.println("Original Count :" +ims1);
//            System.out.println("read add Count :" +iaddCount);

            System.out.println("after combine CheckPeptideWithAllChargeAtSameRT size:"+sorted.size());


           /* if(Utils.bzCalculateThePearsonCorresionWithPredictMSMS)
            {
                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致

                *//*FileReader freader;
                BufferedReader br;
                String dataFile = Config.spectrumFolder +Utils.strpDeepOutfile;//pdeep3output.txt";
                String line;
                try {
                    freader = new FileReader(dataFile);
                    br = new BufferedReader(freader);
                    while ((line = br.readLine()) != null && !(line.startsWith(">peptide|")));//定位到当前位置

                } catch (FileNotFoundException noFile) {
                    throw new FileNotFoundException();
                }*//*


                for(int iS = 0 ;iS<sorted.size();iS++) {
//                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS_JustAll(
                            MapMS2Intensity.get(
                                    sorted.get(iS).strPepcomposition+sorted.get(iS).msOneTrail.getZ()).arrdIntensity);
                    if(iS%10000==0)
                        System.out.println(iS + "---------iS---------" + sorted.get(iS).getdMatch() + sorted.get(iS).strDecoy);

                   *//* //open predict MSMS results file
                    double[][] arrPredictMSMS = new double[sorted.get(iS).pep.composition.length() - 1][4];//length -1 with b b+2 y y+2

                    String[] strPredictPepInfo = line.split("\\|");

                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序

                    while (!strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        if ((line = br.readLine()) != null && (line.startsWith(">peptide|"))) {
                            strPredictPepInfo = line.split("\\|");
                        }
                    }
                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序，过滤掉不匹配的数据

                    if (strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        String strValue = "";
                        while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) {
                            strValue = strValue + line;
                        }
                        String[] arrstrValue = strValue.replaceAll("\\[ ", "")
                                .replaceAll("\\[", "")
                                .replaceAll("\\]", "")
                                .replaceAll("  ", " ")
                                .split(" ");
                        for (int iP = 0; iP < strPredictPepInfo[1].length() - 1; iP++) {
                            for (int jP = 0; jP < 4; jP++) {
                                if(arrstrValue[iP * 4 + jP].isEmpty())
                                {
                                    System.out.println(sorted.get(iS).pep.composition);//test
                                }else
                                {
                                    arrPredictMSMS[iP][jP] = Double.parseDouble(arrstrValue[iP * 4 + jP]);

                                }
                            }
                        }
                    }
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(arrPredictMSMS);*//*
                }
               // br.close();
               // freader.close();
                //predict msms results end
            }
            System.out.println("---------CalculateThePearsonCorresionWithPredictMSMS Finish---------");
*/
            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long,Integer> mapIDCount = new HashMap();

            if (Utils.bzJustOneTrailCanMatchAPeptide) //need to adjust the match score
            {

                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致


                for(int iS = 0 ;iS<20000 && iS<sorted.size();iS++) {
                    int jS = iS + 1;
                    if (sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) {
                                l = mid + 1;
                            } else {
                                h = mid;
                            }
                        }

                        if (l - iS > 1 && sorted.get(iS).getdMatch() < sorted.get(iS + 1).getdMatch()) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS + "----jS--------------" + (jS < sorted.size() ? sorted.get(jS).getdMatch() : "last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        } else {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }
                    } else {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }
                }

                 /*   if(sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) l = mid + 1;
                            else h = mid;
                        }


//                        for(;jS<sorted.size();jS++)
//                        {
//                            if(sorted.get(iS).getdMatch() >= sorted.get(jS).getdMatch()) break;
//                        }
                        if(jS-iS>1) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS+"----jS--------------"+(jS<sorted.size()?sorted.get(jS).getdMatch():"last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        }else
                        {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }

                    }else
                    {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }

                }*/
/*                for (MatchPeptdide mp : sorted) {
                    mp.adjustScore();
                }*/


            }
            Collections.sort(sorted, Collections.reverseOrder());

            //output multiple file, concise file, only peptidefile
            FileWriter concisefileWriter;
            BufferedWriter concisebufferedWriter = null;
            FileWriter peptidesForpDeep3Writer = null;
            BufferedWriter peptidesForpDeep3bufferedWriter = null;
            FileWriter peptidesForpPrositWriter = null;
            BufferedWriter peptidesForpPrositbufferedWriter = null;
            if(Utils.bzOutputConciseAndPeptideFiles)
            {
                concisefileWriter = new FileWriter(Config.spectrumFolder + "concisefile_0502.csv");
                concisebufferedWriter = new BufferedWriter(concisefileWriter);
                peptidesForpDeep3Writer = new FileWriter(Config.spectrumFolder + "peptidesForpDeep3_0502.csv");
                peptidesForpDeep3bufferedWriter = new BufferedWriter(peptidesForpDeep3Writer);
                peptidesForpPrositWriter = new FileWriter(Config.spectrumFolder + "peptidesForpProsit_0502.csv");
                peptidesForpPrositbufferedWriter = new BufferedWriter(peptidesForpPrositWriter);
                concisebufferedWriter.write("msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tgetQuality_score\tmsonemz\tmsonecharge\tBIon\tYIon"+"\n");
                peptidesForpDeep3bufferedWriter.write("peptide\tmodinfo\tcharge"+"\n");
                peptidesForpPrositbufferedWriter.write("modified_sequence\tcollision_energy\tprecursor_charge\tfragmentation\n");

            }
            System.out.println("---------OutputALL Begin---------");
            ProgressBar pboutput = new ProgressBar("OutputALL", sorted.size());
            pboutput.start();
            for (MatchPeptideWithMultiCharge mp : sorted) {
                pboutput.step();
                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
//                bufferedWriter.write(k+ "\t");
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");

//                    if(Utils.bzOutputConciseAndPeptideFiles)
//                    {
//
//                        concisebufferedWriter.write(mp.getConciseString()+"\n");
//                        peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
//                        peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");
//
//                    }
                }else
                {
//                    System.out.println(mp.strPepcomposition);
                    bufferedWriter.write(mp.toString() + "\n");


                }
                if(Utils.bzOutputConciseAndPeptideFiles)
                {

                    concisebufferedWriter.write(mp.getConciseString()+"\n");
                    peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
                    peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");

                }

                iCountForFDR++;

            }
            pboutput.stop();
            System.out.println("---------OutputALL Finish---------");

            if(Utils.bzOutputConciseAndPeptideFiles) {
                concisebufferedWriter.flush();
                peptidesForpDeep3bufferedWriter.flush();
                peptidesForpPrositbufferedWriter.flush();
                concisebufferedWriter.close();
                peptidesForpDeep3bufferedWriter.close();
                peptidesForpPrositbufferedWriter.close();
            }
        }
//        System.out.println("The number of matched peps is: " + );


        bufferedWriter.flush();
        bufferedWriter.close();

        // sort result
//        Collections.sort(resLst, Collections.reverseOrder());
//        System.out.println("The number of matched peps is: " + resLst.size());
//        bufferedWriter.flush();
//        bufferedWriter.close();

//        for (DbMatch match : resLst) {
//            // TODO Remove all System.out.println calls for production
//            System.out.println(match.id);
//            System.out.println(match.composition);
//            System.out.println(match.matchedRts.score);
//
//            // Output the information to a file
//            bufferedWriter.write(match.id + "\n");
//            bufferedWriter.write(match.composition + "\n");
//            bufferedWriter.write(Double.toString(match.matchedRts.score) + "\n");
//
//            if (Config.mode == Enums.RunMode.DEBUG) {
//                System.out.println(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")) + "\n");
//            } else {
//                System.out.println(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")) + "\n");
//            }
//        }
    }

    public void match_TrailWithMultiRT_FitRT_PeptideAllChargeCheckAtsameRT_CorssWindow(String mgfInfile,String rawfile, String fastaInfile, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();

        // read fasta
        ArrayList<Genome> genomes = FastaFile.ReadFile(fastaInfile);
//        ArrayList<Genome> genomes = FastaFile.ReadFileFormPeptide(fastaInfile);


        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top
//        List<Peptide> peps = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
//        double[] pepsMass = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().mapToDouble(x -> x.mass).toArray();

        List<Peptide> pepsOri = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).sorted().collect(Collectors.toCollection(ArrayList::new));


        if (Utils.bRemovePeptideFromBothTargetAndDecoy) {

            TreeMap<String, List<Peptide>> mapPepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(Peptide::getComposition, TreeMap::new, Collectors.toList()));

            TreeMap<String, List<Peptide>> mapILreplacePepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(Peptide::getReplaceILComposition, TreeMap::new, Collectors.toList()));

            System.out.println(mapPepstrListPep.size());
            ArrayList<Peptide> arrPepsRemove = new ArrayList<>();

            for (String strPep : mapPepstrListPep.keySet()) {


                boolean bDecoyType = false;
                boolean bTargetType = false;

                if (mapPepstrListPep.get(strPep).size() > 1) {
                    for (Peptide pep : mapPepstrListPep.get(strPep)) {

                        boolean isDecoy =pep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }
                }
                //REMOVE I L MUTATE peptide in both decoy and target
                if(strPep.contains("I")||strPep.contains("L"))
                {
                    for(Peptide muSamePep:mapILreplacePepstrListPep.get(strPep.replaceAll("L","I")))
                    {
                        boolean isDecoy =muSamePep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }

                }


                if (!(bDecoyType && bTargetType)) {
                    arrPepsRemove.addAll(mapPepstrListPep.get(strPep));
                }

            }
            System.out.println("The number of arrPepsRemove is: " + arrPepsRemove.size());

            pepsOri = arrPepsRemove;

        }
        List<Peptide> peps = pepsOri.stream().distinct().sorted().collect(Collectors.toCollection(ArrayList::new));



        double[] pepsMass = peps.stream().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());
        pepsOri = null;


//        FileWriter fileWriterPeptideInfo = new FileWriter("peptideInfo.csv");
//        BufferedWriter bufferedWriterPeptideInfo = new BufferedWriter(fileWriterPeptideInfo);
//        bufferedWriterPeptideInfo.write( "peptide\tmass\tdMutationRate"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
//                + "\n");
//        for (Peptide pep:peps)
//        {
//            bufferedWriterPeptideInfo.write(pep.composition+'\t'+pep.mass+'\t'+pep.dMutationRate+'\n');
//        }
//
//        bufferedWriterPeptideInfo.flush();
//        bufferedWriterPeptideInfo.close();


        ArrayList<DbMatch> resLst = new ArrayList<>();




//        Set<Peptide>
        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tmsoneMZ\tmsoneZ\tmsonePeakAreaLocalRank" +
                "\tms2rt\tpeptideLength\tMS1PeptideMassError\tdPeptideMutationRate\tibyComplemetaryCount\tdCosWithPredictMSMS\tdMS1WithMS2TrailCoscinSimilarity\tMS2TrailCoscinSimilarity\tBIon\tYIon");
        for(int iw =0;iw<Utils.iMS2Charge;iw++) {
            bufferedWriter.write("\tdlogBionPeakAreaSum_C"+(iw+1)+"\tdlogYionPeakAreaSum_C"+(iw+1)+"\tdBionMassErrorSum_C"+(iw+1)+"\tdYionMassErrorSum_C"+(iw+1)+"\tiBionPeakHalfSum_C"+(iw+1)+"\tiYionPeakHalfSum_C"+(iw+1)+"\tbionMatchCount_C"+(iw+1)+"\tyionMatchCount" +
                    "_C"+(iw+1)+"\tdBionRetentionTimeErrorSum_C"+(iw+1)+"\tdYionRetentionTimeErrorSum_C"+(iw+1)+"\tiBConsective_C"+(iw+1)+"\tiYConsective_C"+(iw+1)+"\tdBionCosWithPredictMSMS_C"+(iw+1)+"\tdYionCosWithPredictMSMS_C"+(iw+1));
        }
        bufferedWriter.write( "\tarrMatchedBion\tarrMatchedYion\tadjustBYIon\tSoredIndex\twindowSize\tdBionCos60WithPredictMSMS_C1\tdCosWithPredictMSMSMatchedPredict\tprositRT" +
                "\tibAllMatchedWithloss\tibMatchNormalAndH2Oloss\tibMatchNormalAndNH3loss\tiyAllMatchedWithloss\tiyMatchNormalAndH2Oloss\tiyMatchNormalAndNH3loss\tisRawdata\tpeptideAnomiAcidFreRate\tlastIonNumBC1\tlastIonNumYC1\tcount"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
                + "\n");



//        Map<Integer,Map<Integer,Double>> mapTriplePeptideRTScore = new TreeMap<>();

//        BitSet mapTriplePeptideRTScore = new BitSet(peps.size());
//        BitSet[] mapTriplePeptideRTScore = new BitSet[peps.size()];
//        Map<Pair<Integer,Integer>,Double> mapTriplePeptideRTScore = new HashMap<>();
//        for (MSTwoTrail msTwoTrail : MSTWOts.arrMSTwoTrail) {
//            double[] rts = msTwoTrail.getRts();
//
//            for(double drt:rts) {
//                if (uniqueValues.add(drt))
//                {
//                    List<MSTwoTrail> listMSTwoTrail = new ArrayList<>();
//                    listMSTwoTrail.add(msTwoTrail);
//
//                    mapRTMSTwoTrail.put(drt,listMSTwoTrail);
//                }
//                else
//                {
//                    mapRTMSTwoTrail.get(drt).add(msTwoTrail);
//                }
//            }
//        }

        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if(!Utils.bzOnlyRawdata)
            {
                spec.readFromeTrailFile(file,false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile,true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


//            icount = spec.arrMSOneTrail.size();
            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();

//            Map<Double, List<MSOneTrail>> mapRTMSoneTrail =  spec.arrMSOneTrail.stream().
//                    collect(Collectors.groupingBy(x-> x.getRt()));

            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));



            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if(RTWidnowDistance!=null && RTWidnowDistance.length>0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
//            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection();
            iswc.IsolationWindowCollection_paralle(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);

//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2Intensity();
//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2IntensityFromProsit();
//            Map<String, MS2Intensity> MapMS2Intensity = new HashMap<>();//Main.getMapMS2IntensityFromProsit();//firstrunn 0829
            Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit();//firstrunn 0829

            //GET The autort information
//            Map<String, Double> MapPeptidePredictAutoRT = Main.getAutoRTInfo();
//            List<String> strpepsAutoRT = MapPeptidePredictAutoRT.keySet().stream().collect(Collectors.toList());
//            double[] dpepsAutoRT = MapPeptidePredictAutoRT.entrySet().stream().mapToDouble(x -> x.getValue()).toArray();
//
//            MapPeptidePredictAutoRT = null; //release



            //set ms1feature peakarea rank
            AtomicInteger ims1 = new AtomicInteger();


            ProgressBar pbMS1 = new ProgressBar("Progress MS1 area rank", arrRts.length);
            pbMS1.start();
//            for(MSOneTrail ms1Trail:spec.arrMSOneTrail)
            spec.arrMSOneTrail.stream().parallel().forEach(ms1->
            {
                MSOneTrail ms1Trail = ms1;
                pbMS1.step();
                long iPeakAreaRank = 1;//至少排名第一
                List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(ms1Trail.getMz(), 1,true);
                for(MSTwoTrailSet  window:arrwindow) {
                    //get ms1feature RT distance between Math.abs(y-v) <= Utils.thresholdRT
                    int iMS1RTPosBegin = Arrays.binarySearch(arrRts, ms1Trail.getRt() - Utils.thresholdMS1CheckPeakAreaRT);//正负15秒内
                    iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
//                if (iMS1RTPosBegin > 0)
//                    iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position

                    for (int iRTPos = iMS1RTPosBegin;
                         iRTPos < arrRts.length && arrRts[iRTPos] < ms1Trail.getRt() + Utils.thresholdMS1CheckPeakAreaRT;
                         iRTPos++) {
                        if (window != null) {
                            //0806加入电荷相等
//                            iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
//                                    x.getZ() == ms1Trail.getZ() && x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh &&
//                                            x.getQuantification_peaks_area() > ms1Trail.getQuantification_peaks_area()
//                            ).count();
                            //0914 multi thread
                            iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
                                    x.getZ() == ms1Trail.getZ() && x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh &&
                                            x.getQuantification_peaks_area() > ms1Trail.getQuantification_peaks_area()
                            ).count();
                        }
                    }
                }
                ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));
//                if(iPeakAreaRank>0)
//                {
//                    ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));
//                }else
//                {
//                    ms1Trail.setdPeakAreaLocalRankInTimeSpan( Math.log(20.0));
//                }


/*                if(window!=null) {
                    iPeakAreaRank = spec.arrMSOneTrail.stream().filter((x) ->
                            x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh && x.getQuantification_peaks_area() >= ms1Trail.getQuantification_peaks_area()
                                    && x.getRt() >= ms1Trail.getRt() - Utils.thresholdMS1RT && x.getRt() <= ms1Trail.getRt() + Utils.thresholdMS1RT
                    ).count();
                }*/

//                System.out.println("ims1:"+ ims1++  +"----"+spec.arrMSOneTrail.size());



            });
            pbMS1.stop();

//        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(
//                "isolationWindowRanges.txt"));
//        IsolationWindowCollection iswc = (IsolationWindowCollection) ois.readObject();

//        List<double[]> larrMSTwoRts = new ArrayList<>(iswTest.windows.size());
//        List<Map<Double, List<MSTwoTrail>>> lmapRTMStwoTrail = new ArrayList<>(iswTest.windows.size());

            // TODO: READ IN SPECTRUM RANKING FILE (Added by Caroline) Implement Rank Here


            // for each isolation window collection, match once
            //Step 1: get a msonetrail
            //Step 2: sort with RT
            //Step 3:for each msonetrail get the window
            //      check all peptide with same mass  (mz +- 30 ppm)
            //         for each mstwoTrail in window (mz in window && RT +- 5s )
            //             compare mstwotrail b y ions with peptide b y ions (mz +- 30 ppm)
            //         add peptide to msonetrail as candidate peptide with confident score
//        MSTwoTrailSet MSTWOts = new MSTwoTrailSet();
//        MSTWOts.readFromeTrailFile(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window0_ms2_trails.tsv");

//        Map<Double,List<MSTwoTrail> >mapRTMSTwoTrail = new HashMap<Double, List<MSTwoTrail>>();
            TreeSet<Double> uniqueValues = new TreeSet<>();


            //////////////////////////////////
            //for a window,initial rts and rts corresponding trails
//        double[] arrMSTwoRts = MSTWOts.arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
//                .distinct()
//                .sorted()
//                .toArray();
////        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().
////                collect(Collectors.groupingBy(x-> x.getRtApex()));
//
//
//        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
            //////////////////////////////////////

//            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = new HashMap<>();
//            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = Collections.synchronizedMap( new ConcurrentHashMap<>());
            List<MatchPeptideWithMultiCharge> listPeptideAdd = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());
            List<MatchPeptideWithMultiCharge> listPeptideAddTrue = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());


//            Map<Long, LinkedList<MatchPeptdide>> mapMSOneFeaturePeptide = new HashMap<>();

//            AtomicInteger iPrecursorNO= new AtomicInteger(1000000);
//            AtomicInteger i = new AtomicInteger();
            ProgressBar pb = new ProgressBar("Progress", arrRts.length);
            pb.start();
//            List<Double> dRTSynlist = Collections.synchronizedList(new ArrayList<>());
//            for(double dRT:arrRts)
//                dRTSynlist.add(dRT);

//            dRTSynlist.stream().parallel().forEach(
//            Arrays.stream(arrRts).forEach(
            Arrays.stream(arrRts).parallel().forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();
                        pb.step();
                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);

                        Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());


                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {
//                                        if(msonetrail.bAddPrecursor)
//                                        {


//                                        if (msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass >= pepsMass[0]) //filter m1 trail which less than the least peptide mass

                                        //begin of window with mz of m1 feature
////                                        MSTwoTrailSet window = iswc.FindWindow(msonetrail.getMass(), 2);
//                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                        // if mass matches, then pep match
//                                        if (window != null) {
//                                            double[] arrMSTwoRts = window.arrMSTwoRts;
//                                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
                                        if (msonetrail.getMass() - Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {
                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window

                                            LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide = new LinkedList<MatchPeptideWithMultiCharge>();

//                                            LinkedList<MatchPeptdide> listMaxiumPeptide = new LinkedList<MatchPeptdide>();


                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况
//                                            while (iPosBegin > 0) {
//                                                iPosBegin = iPosBegin - 1;//move to less than position
//                                                if(iPosBegin > 0 && pepsMass[iPosBegin]>pepsMass[iPosBegin-1]) break;//PEPTIDE的质量可能有相等的情况，这是peptide有序情况
//                                            }

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {



                                                //begin of window with mz of ms1 feature
//                                                MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
                                                // if mass matches, then pep match
                                                List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(msonetrail.getMz(), 1,true);
                                                for(MSTwoTrailSet  window:arrwindow) {
                                                    MS1MatchPeptide(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msonetrail, listMaxiumPeptide, peps.get(iPos), window);
                                                }
                                            }

//                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);

                                            if (listMaxiumPeptide.size() > 0)
                                            {
                                                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                                                {
                                                    listPeptideAddTrue.addAll(listMaxiumPeptide);

                                                }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
                                                {
                                                    listPeptideAddTrue.add(listMaxiumPeptide.get(0));
                                                }
                                                else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
                                                {
                                                    if (listMaxiumPeptide.size() > 1) {
                                                        listPeptideAddTrue.add(listMaxiumPeptide.get(1));
                                                    }
                                                }
                                            }
//                                            listPeptideAddTrue.addAll(listMaxiumPeptide);
                                            ims1.addAndGet(listMaxiumPeptide.size());
                                            if (Utils.bzCheckPeptideWithAllChargeAtSameRT)
                                            {
                                                //检查该MS1相同RT且已经匹配的相同PEPTIDE里面所有电荷的二级谱
                                                int iCurrenMS1Z = msonetrail.getZ();
//                                                Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());
                                                for(MatchPeptideWithMultiCharge matchPep:listMaxiumPeptide)
                                                {
                                                    //是否重复出现的PEPTIDE，若相同只进来一次
                                                    if(uniquePeptideString.contains(matchPep.pep.composition)) {
                                                        continue;
                                                    }
                                                    for(int iAllC=1;iAllC<=Utils.iMS1Charget;iAllC++)
                                                    {
                                                        //condition is the z is larger than current z
                                                        if(iAllC == iCurrenMS1Z) {
                                                            continue;//means this
                                                        }
                                                        double dPepMZWithC = Utils.MassToMz(matchPep.pep.mass,iAllC);
                                                        //search the mz in lMsoneTrail,如果存在（可能存在一个或多个）且比当前小就不search 二级谱，否则search，然后在标记该peptide已经被search，下次就不再search

                                                        if(iAllC<iCurrenMS1Z) {
                                                            //if it is matched, it should told the problem there is no matched at successive process(create a map??)
                                                            //get the window
                                                            boolean isE = false;
                                                            int iMS1PosBegin = Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            for (int iPos = iMS1PosBegin;
                                                                 iPos < lMsoneTrail.size() && lMsoneTrail.get(iPos).getMz() < dPepMZWithC*(1+Utils.thresholdPPM);//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                                 iPos++) {

                                                                isE =true;
                                                                break;

                                                            }
                                                            if(isE){
                                                                continue;
                                                            }

                                                        }else if (iAllC>iCurrenMS1Z){
                                                            //(the larger z means mz is less than current mz that has been checked before)
                                                            int iMS1PosBegin = Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC*(1-Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            if(iMS1PosBegin
                                                                    < lMsoneTrail.size() && lMsoneTrail.get(iMS1PosBegin).getMz() < dPepMZWithC*(1+Utils.thresholdPPM))//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                            {

                                                                continue;
                                                            }


                                                        }
//                                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(dPepMZWithC, 1);
                                                        List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(dPepMZWithC, 1,true);
                                                        for(MSTwoTrailSet  window:arrwindow) {
                                                            MSOneTrail msOneTrail = new MSOneTrail(-1, dPepMZWithC, msonetrail.getRt(), iAllC, 0, 0, 0);
                                                            LinkedList<MatchPeptideWithMultiCharge> listPeptide = new LinkedList<MatchPeptideWithMultiCharge>();
                                                            MS1MatchPeptide(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msOneTrail, listPeptide, matchPep.pep, window);
                                                            listPeptideAdd.addAll(listPeptide);
                                                        }

                                                        //check other charge



                                                    }
                                                    //check whether it has matched before

                                                    uniquePeptideString.add(matchPep.pep.composition);

                                                }


                                            }
                                        }
//                                        }


                                        //get all peps mass equals the msonetrail
//                                        peps.stream()
//                                                .filter(e->Math.abs(e.mass-msonetrail.getMass())<Utils.thresholdPPM*e.mass)
//                                                .forEach(
//                                                        pep->{
//                                                            //compare pep with all mstwotrail in rt window by using b y ions
//                                                            uniqueValues.stream().filter(v->Math.abs(y-v) <= Utils.thresholdRT).forEach(
//                                                                    twoRTtime->{
////                                                                         if(DbMatch.PepMatch(pep,mapRTMSTwoTrail.get(twoRTtime),Config.scoreMode)>0.0) {
////                                                                             System.out.println(icount.incrementAndGet() +"--------");
////                                                                             System.out.println(pep.composition + pep.mass);
////                                                                             System.out.println(msonetrail.getMass());
////                                                                             System.out.println(twoRTtime);
////                                                                         }
//
////                                                                        mapRTMSTwoTrail.get(twoRTtime).stream().forEach(
////                                                                                msTwoTrail -> {
////                                                                                    DbMatch.PepMatch(pep,msTwoTrail,Config.scoreMode);
////
////
////                                                                                }
////                                                                        );
//                                                                    }
//                                                            );
//
//
//                                                        }
//                                                );
//                                    }//end of the window with mz of ms1 feature
                                    }
                            );
                        }





//                        i.getAndIncrement();
//                        long time = System.currentTimeMillis() - start;
                        //                       System.out.println("#########time:" + i + "/" + arrRts.length);

                        //                       System.out.println(time);
                    }
            );
            pb.stop();




//            //根据RT，和预测的RT的拟合值，找出在RT范围内的PEPTIDE，添加precursor到现有的msoneTRAIL中。
//            //设置这个precursor只检查这个PEPTIDE(以前是检查所有peptide，导致非常慢，可以将rt慢长一些[-6.5,+6.5]，加1-4个电荷)
//            //考虑不同的电荷检查raw data里面是否存在误差容限内的isotope是否存在，可以检查一到两个isotope：1电荷是+1，2电荷是+0.5，3电荷是+0.33，4电荷是+0.25
//            //get  RT distance between ms1 interval 得到一个rt间隔的所有PEPTIDE
//
////            int iAutoRTPosBegin = Arrays.binarySearch(dpepsAutoRT, y-RTWidnowDistance[0]);
////            iAutoRTPosBegin = Utils.getiPosBeginOfNearest(dpepsAutoRT.length, iAutoRTPosBegin);
////            if (iAutoRTPosBegin > 0)
////                iAutoRTPosBegin = iAutoRTPosBegin - 1;//move to less than position
//            for (int iRTPos = 0;
//                 iRTPos < dpepsAutoRT.length ;
//                 iRTPos++) {
//
//                String strPep = strpepsAutoRT.get(iRTPos);
//                double dmass = Peptide.CalculateMassByPeptideStr(strPep);
//                for (int icharge = 1;icharge<=Utils.iChargeRange;icharge++)
//                {
//                    double dmz = Utils.MassToMz(dmass,icharge);
//
//                    //直接匹配相对应的PEPTIDE
//                    MSTwoTrailSet window = iswc.FindWindowWithMZ(dmz, 1);
////                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
//                    // if mass matches, then pep match
//                    if (window != null) {
//                        double[] arrMSTwoRts = window.arrMSTwoRts;
//                        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
//
//                        double rtbegin = dpepsAutoRT[iRTPos]-Utils.thresholdPredictRT;
//
//                        int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, rtbegin);
//                        iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
//                        if (iTwoRTPosBegin > 0)
//                            iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
//
//
//                        for (int iMS2RTPos = iTwoRTPosBegin;
//                             iMS2RTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < dpepsAutoRT[iRTPos] + Utils.thresholdPredictRT;
//                             iMS2RTPos++) {
//                            //只匹配一个PEPTIDE,但基于多个时间
//                            peps.get(iPos).GenerateIons();
//                            double[][] arrIntensity = null;
//                            MS2Intensity ms2Int = MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ());
//                            if (ms2Int!=null)
//                                arrIntensity = ms2Int.arrdIntensity;
//
//                            MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiCharge(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                    Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail
//                                    ,arrIntensity);
//                                                       /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
//                                                                msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
//                            if (matchPep != null) {
//                                matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错
////                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);
//
//                                if (ms2Int != null)
//                                    matchPep.dPrositRT = ms2Int.dPredictRT;
//                                matchPep.calculateCombineScore();
//                                matchPep.setDms2rt(arrMSTwoRts[iRTPos]);
//
//
//                                matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);
//                            }
//
//                        }
////                        MSOneTrail tmpPrecursor = new MSOneTrail(
////                                iPrecursorNO.getAndIncrement(), dmz, dpepsAutoRT[iRTPos],
////                                icharge, 0, 0,
////                                0, 0, 0,
////                                0, 0, 0,
////                                0, String.valueOf(dmz), String.valueOf(dpepsAutoRT[iRTPos]), "0"
////
////                        );
////                        tmpPrecursor.strGeneratePeptide = strPep;
////                        tmpPrecursor.bAddPrecursor = true;//set the precurosor bz
////                        lMsoneTrail.add(tmpPrecursor);
//                    }
//                }
//
//            }

//            for (Peptide pep: peps) {
//                // search for the window
//                IsolationWindow window = spec.FindWindow(pep.mass, 2);
//                // if mass matches, then pep match
//                if (window != null) {
//                    DbMatch res = DbMatch.PepMatch(pep, window, Config.scoreMode);
//                    if (res != null) {
//                        if (res.matchedRts.score > 0) {
//                            resLst.add(res);
//                        }
//                    }
//                }
//            }




//            List<MatchPeptdide> sorted = new LinkedList<>();
            System.out.println("add listPeptideAddTrue size:"+listPeptideAddTrue.size());

            System.out.println("add CheckPeptideWithAllChargeAtSameRT size:"+listPeptideAdd.size());
//            List<MatchPeptideWithMultiCharge> sorted = new LinkedList<>();
            List<MatchPeptideWithMultiCharge> sorted = listPeptideAdd;
            sorted.addAll(listPeptideAddTrue);
//            for(MatchPeptideWithMultiCharge map:sorted)
//            {
//                if (map==null)
//                {
//                    System.out.println("null");
//                }
//            }


//            List<MatchPeptdide> sorted = mapMSOneFeaturePeptide.entrySet().stream()
////                    .sorted(Comparator.comparing(e -> e.getValue().stream().map(MatchPeptdide::getdMatch).min(Comparator.reverseOrder()).orElse((double) 0)))
//                    //and also sort each group before collecting them in one list
//                    .flatMap(e -> e.getValue().stream().sorted(Comparator.comparing(MatchPeptdide::getdMatch))).collect(Collectors.toList());
//
//            for (MatchPeptdide mp:sorted)
//                {
//                    bufferedWriter.write(mp.toString()+"\n");
//
//                }
            //
            //
//            Set<Long> keys = mapMSOneFeaturePeptide.keySet();

//            System.out.println("matched MSonePrecursor Size:"+keys.size());
//            if (Utils.bzJustOneTrailCanMatchAPeptide) Utils.OutputMode = Utils.OutputModeEnum.OutputALL;//consider all ms1feature can match peptide, because those match score should be adjust

//            int iaddCount=0;
//            int iPrecursorCount = 0;
//            for (Long k : keys) {
//
//                List<MatchPeptideWithMultiCharge> lmp = mapMSOneFeaturePeptide.get(k);
////                List<MatchPeptdide> lmp = mapMSOneFeaturePeptide.get(k);
//                iPrecursorCount++;
//                if (lmp.size() > 0)
//                {
//                    if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
//                    {
//                        sorted.addAll(lmp);
//
//                        iaddCount+=lmp.size();
//
//                    }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
//                    {
//                        sorted.add(lmp.get(0));
//                    }
//                    else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
//                    {
//                        if (lmp.size() > 1) {
//                            sorted.add(lmp.get(1));
//                        }
//                    }
//                }
//
////                sorted.addAll(lmp);
//
//            }
//            System.out.println("precursor Count :" +iPrecursorCount);

            System.out.println("Original Count :" +ims1);
//            System.out.println("read add Count :" +iaddCount);

            System.out.println("after combine CheckPeptideWithAllChargeAtSameRT size:"+sorted.size());


           /* if(Utils.bzCalculateThePearsonCorresionWithPredictMSMS)
            {
                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致

                *//*FileReader freader;
                BufferedReader br;
                String dataFile = Config.spectrumFolder +Utils.strpDeepOutfile;//pdeep3output.txt";
                String line;
                try {
                    freader = new FileReader(dataFile);
                    br = new BufferedReader(freader);
                    while ((line = br.readLine()) != null && !(line.startsWith(">peptide|")));//定位到当前位置

                } catch (FileNotFoundException noFile) {
                    throw new FileNotFoundException();
                }*//*


                for(int iS = 0 ;iS<sorted.size();iS++) {
//                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS_JustAll(
                            MapMS2Intensity.get(
                                    sorted.get(iS).strPepcomposition+sorted.get(iS).msOneTrail.getZ()).arrdIntensity);
                    if(iS%10000==0)
                        System.out.println(iS + "---------iS---------" + sorted.get(iS).getdMatch() + sorted.get(iS).strDecoy);

                   *//* //open predict MSMS results file
                    double[][] arrPredictMSMS = new double[sorted.get(iS).pep.composition.length() - 1][4];//length -1 with b b+2 y y+2

                    String[] strPredictPepInfo = line.split("\\|");

                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序

                    while (!strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        if ((line = br.readLine()) != null && (line.startsWith(">peptide|"))) {
                            strPredictPepInfo = line.split("\\|");
                        }
                    }
                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序，过滤掉不匹配的数据

                    if (strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        String strValue = "";
                        while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) {
                            strValue = strValue + line;
                        }
                        String[] arrstrValue = strValue.replaceAll("\\[ ", "")
                                .replaceAll("\\[", "")
                                .replaceAll("\\]", "")
                                .replaceAll("  ", " ")
                                .split(" ");
                        for (int iP = 0; iP < strPredictPepInfo[1].length() - 1; iP++) {
                            for (int jP = 0; jP < 4; jP++) {
                                if(arrstrValue[iP * 4 + jP].isEmpty())
                                {
                                    System.out.println(sorted.get(iS).pep.composition);//test
                                }else
                                {
                                    arrPredictMSMS[iP][jP] = Double.parseDouble(arrstrValue[iP * 4 + jP]);

                                }
                            }
                        }
                    }
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(arrPredictMSMS);*//*
                }
               // br.close();
               // freader.close();
                //predict msms results end
            }
            System.out.println("---------CalculateThePearsonCorresionWithPredictMSMS Finish---------");
*/
            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long,Integer> mapIDCount = new HashMap();

            if (Utils.bzJustOneTrailCanMatchAPeptide) //need to adjust the match score
            {

                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致


                for(int iS = 0 ;iS<20000 && iS<sorted.size();iS++) {
                    int jS = iS + 1;
                    if (sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) {
                                l = mid + 1;
                            } else {
                                h = mid;
                            }
                        }

                        if (l - iS > 1 && sorted.get(iS).getdMatch() < sorted.get(iS + 1).getdMatch()) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS + "----jS--------------" + (jS < sorted.size() ? sorted.get(jS).getdMatch() : "last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        } else {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }
                    } else {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }
                }

                 /*   if(sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) l = mid + 1;
                            else h = mid;
                        }


//                        for(;jS<sorted.size();jS++)
//                        {
//                            if(sorted.get(iS).getdMatch() >= sorted.get(jS).getdMatch()) break;
//                        }
                        if(jS-iS>1) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS+"----jS--------------"+(jS<sorted.size()?sorted.get(jS).getdMatch():"last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        }else
                        {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }

                    }else
                    {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }

                }*/
/*                for (MatchPeptdide mp : sorted) {
                    mp.adjustScore();
                }*/


            }
            Collections.sort(sorted, Collections.reverseOrder());

            //output multiple file, concise file, only peptidefile
            FileWriter concisefileWriter;
            BufferedWriter concisebufferedWriter = null;
            FileWriter peptidesForpDeep3Writer = null;
            BufferedWriter peptidesForpDeep3bufferedWriter = null;
            FileWriter peptidesForpPrositWriter = null;
            BufferedWriter peptidesForpPrositbufferedWriter = null;
            if(Utils.bzOutputConciseAndPeptideFiles)
            {
                concisefileWriter = new FileWriter(Config.spectrumFolder + "concisefile_0502.csv");
                concisebufferedWriter = new BufferedWriter(concisefileWriter);
                peptidesForpDeep3Writer = new FileWriter(Config.spectrumFolder + "peptidesForpDeep3_0502.csv");
                peptidesForpDeep3bufferedWriter = new BufferedWriter(peptidesForpDeep3Writer);
                peptidesForpPrositWriter = new FileWriter(Config.spectrumFolder + "peptidesForpProsit_0502.csv");
                peptidesForpPrositbufferedWriter = new BufferedWriter(peptidesForpPrositWriter);
                concisebufferedWriter.write("msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tgetQuality_score\tmsonemz\tmsonecharge\tBIon\tYIon"+"\n");
                peptidesForpDeep3bufferedWriter.write("peptide\tmodinfo\tcharge"+"\n");
                peptidesForpPrositbufferedWriter.write("modified_sequence\tcollision_energy\tprecursor_charge\tfragmentation\n");

            }
            System.out.println("---------OutputALL Begin---------");
            ProgressBar pboutput = new ProgressBar("OutputALL", sorted.size());
            pboutput.start();
            for (MatchPeptideWithMultiCharge mp : sorted) {
                pboutput.step();
                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
//                bufferedWriter.write(k+ "\t");
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");

//                    if(Utils.bzOutputConciseAndPeptideFiles)
//                    {
//
//                        concisebufferedWriter.write(mp.getConciseString()+"\n");
//                        peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
//                        peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");
//
//                    }
                }else
                {
//                    System.out.println(mp.strPepcomposition);
                    bufferedWriter.write(mp.toString() + "\n");


                }
                if(Utils.bzOutputConciseAndPeptideFiles)
                {

                    concisebufferedWriter.write(mp.getConciseString()+"\n");
                    peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
                    peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");

                }

                iCountForFDR++;

            }
            pboutput.stop();
            System.out.println("---------OutputALL Finish---------");

            if(Utils.bzOutputConciseAndPeptideFiles) {
                concisebufferedWriter.flush();
                peptidesForpDeep3bufferedWriter.flush();
                peptidesForpPrositbufferedWriter.flush();
                concisebufferedWriter.close();
                peptidesForpDeep3bufferedWriter.close();
                peptidesForpPrositbufferedWriter.close();
            }
        }
//        System.out.println("The number of matched peps is: " + );


        bufferedWriter.flush();
        bufferedWriter.close();

        // sort result
//        Collections.sort(resLst, Collections.reverseOrder());
//        System.out.println("The number of matched peps is: " + resLst.size());
//        bufferedWriter.flush();
//        bufferedWriter.close();

//        for (DbMatch match : resLst) {
//            // TODO Remove all System.out.println calls for production
//            System.out.println(match.id);
//            System.out.println(match.composition);
//            System.out.println(match.matchedRts.score);
//
//            // Output the information to a file
//            bufferedWriter.write(match.id + "\n");
//            bufferedWriter.write(match.composition + "\n");
//            bufferedWriter.write(Double.toString(match.matchedRts.score) + "\n");
//
//            if (Config.mode == Enums.RunMode.DEBUG) {
//                System.out.println(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")) + "\n");
//            } else {
//                System.out.println(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")) + "\n");
//            }
//        }
    }
    public void match_FDRLT01_CheckSharePeaksAtRT1MinRangeMT_CorssWindow(String mgfInfile,String rawfile, String fastaInfile,String resultPSMFileName, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();



        //read file from top1 fdr<0.01 dataset
        ArrayList<PSMResult> arrPsmResults = SharePeaks.ReadPSMFromResultFile(resultPSMFileName);


        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);



        bufferedWriter.write("peptide\tmsonefeature\tproteinname\tdecoyTarget\tpepmass\tmsonemass\tmsonetime\tmsoneMZ" +
                "\tmsoneZ\tBIon\tYIon\ttargetscore\tcountDP\tfdr\tsharepeaksMAX\tsharePeaksAll" +
                "\tstrSharePeakMax\tstrSharePeakAll\tshareMs1feature\tshareMS1MZ\tshareMS1time\n");



        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if (!Utils.bzOnlyRawdata) {
                spec.readFromeTrailFile(file, false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile, true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();


            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));


            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if (RTWidnowDistance != null && RTWidnowDistance.length > 0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
//            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName, mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection();
            iswc.IsolationWindowCollection_paralle(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);


            ProgressBar pb = new ProgressBar("Progress", arrPsmResults.size());
            pb.start();

            Map<PSMResult,MatchPeptideWithMultiCharge> mapPSMBYIon = Collections.synchronizedMap( new ConcurrentHashMap<>());//new HashMap<>();


            //the first psm
//            arrPsmResults.get(0).sharepeaksMax = 0;
//            MSTwoTrailSet windowTop = iswc.FindWindowWithMZ( arrPsmResults.get(0).msoneMZ, 1);
//            double[] arrMSTwoRts = windowTop.arrMSTwoRts;
//            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = windowTop.mapRTMStwoTrail;
//            MatchPeptideWithMultiCharge matchHighPSM = getBYIonsMatchWithPSM(arrPsmResults.get(0), windowTop, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
//            matchHighPSM.setMS2TrailSelected();
//            mapPSMBYIon.put(arrPsmResults.get(0),matchHighPSM);
//            bufferedWriter.write(arrPsmResults.get(0) +"\t0\t0\t0\n");

            arrPsmResults.stream().parallel().forEach(
                    y -> {
                        pb.step();
//                        MSTwoTrailSet window = iswc.FindWindowWithMZ( y.msoneMZ, 1);
                        double dmatch = -100000000;
                        List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(y.msoneMZ, 1,true);
//                        if(y.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
//                        {
//                            System.out.println(arrwindow.size());
//                        }
                        for(MSTwoTrailSet  window:arrwindow) {
                            double[] arrMSTwoRts = window.arrMSTwoRts;
                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
                            MatchPeptideWithMultiCharge matchPSM = getBYIonsMatchWithPSM(y, window, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail, arrRts);
//                            if (matchPSM == null)
//                                System.out.println(y.strpeptide + " " + y.msoneZ + " " + y.msonetime + " " + y.msoneMZ);
                            if(matchPSM != null) {
//                                if(y.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
//                                {
//                                    System.out.println(matchPSM.dMatch);
//                                }
                                if(matchPSM.dMatch>dmatch) {//save the max MATCH
                                    mapPSMBYIon.put(y, matchPSM);
                                    dmatch = matchPSM.dMatch;
                                }
                            }
                        }
                    }
            );
            pb.stop();

            List<PSMResult> arrPsm = mapPSMBYIon.keySet().stream().sorted(Comparator.comparingDouble(PSMResult::getTargetscore).reversed()).collect(Collectors.toList());
            ProgressBar pbafter = new ProgressBar("Progress after", arrPsm.size());
            pbafter.start();

            for (int i = 0; i < arrPsm.size(); i++) {
                pbafter.step();
//                if(i%100000==0)
//                    System.out.println(i+"/"+arrPsm.size());
                //get rt range[-1min,1min] and the same mz window from [0,i-1]
                PSMResult curPsm = arrPsmResults.get(i);

                int iSharePeak = 0;
                double maxShareMS1MZ = 0.0;
                double maxShareMS1Time = 0.0;
                String strshareMax="";
                Long ims1id =0L;
                for (int is = 0; is < i; is++) {

                    PSMResult highPsm = arrPsmResults.get(is);
                    if (curPsm.msonetime > highPsm.msonetime - 1 && curPsm.msonetime < highPsm.msonetime + 1
                            && iswc.isSmaeMS1Window(curPsm.msoneMZ, highPsm.msoneMZ)) {

//                        MSTwoTrailSet window = iswc.FindWindowWithMZ(curPsm.msoneMZ, 1);

//                        int iCurSharePeak = getSharePeak(curPsm, highPsm,window,mapRTMSoneTrail,mapPSMBYIon,arrRts);
                        int iCurSharePeak = getSharePeak(curPsm, highPsm,mapPSMBYIon);
                        if (iCurSharePeak > iSharePeak)
                        {
                            iSharePeak = iCurSharePeak;
                            maxShareMS1MZ = highPsm.msoneMZ;
                            maxShareMS1Time = highPsm.msonetime;
                            ims1id = highPsm.msonefeature;
                            strshareMax = mapPSMBYIon.get(curPsm).strSharePeaksMaxInfo;
                        }
                    }


                }
                curPsm.sharepeaksMax = iSharePeak;
                curPsm.strSharePeaksMaxInfo = strshareMax;
                if(mapPSMBYIon.get(curPsm)!=null) {
//                    curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALL();
                    curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALLSpeed();
                    curPsm.strSharePeaksAllInfo = mapPSMBYIon.get(curPsm).strSharePeaksAllInfo;
                    mapPSMBYIon.get(curPsm).setMS2TrailSelectedSpeed();
//                    mapPSMBYIon.get(curPsm).setMS2TrailSelected();


                }
//                else {
//                    MSTwoTrailSet windowCurPSM = iswc.FindWindowWithMZ( curPsm.msoneMZ, 1);
//                    double[] arrMSTwoRtsCurPSM = windowCurPSM.arrMSTwoRts;
//                    Map<Double, List<MSTwoTrail>> mapRTMStwoTrailCurPSM = windowCurPSM.mapRTMStwoTrail;
//                    MatchPeptideWithMultiCharge matchHighPSMCurPSM = getBYIonsMatchWithPSM(curPsm, windowCurPSM, mapRTMSoneTrail, arrMSTwoRtsCurPSM, mapRTMStwoTrailCurPSM,arrRts);
//                    matchHighPSMCurPSM.setMS2TrailSelected();
//                    mapPSMBYIon.put(curPsm,matchHighPSMCurPSM);
//                }

                bufferedWriter.write(curPsm +"\t"+ims1id+"\t"+maxShareMS1MZ+"\t"+maxShareMS1Time+"\n");


            }
            pbafter.stop();


        }


        bufferedWriter.flush();
        bufferedWriter.close();

    }

    public void match_FDRLT01_CheckSharePeaksAtRT1MinRangeMT_CorssWindow_parallel_withMultiplePSMFile
            (String mgfInfile,String rawfile, String fastaInfile,String spectrumFolder ,String resultPSMFileName, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        String[] psmfiles = resultPSMFileName.split(",");
        String[] psmfilesOuts = psmOutfile.split(",");

        AtomicInteger icount = new AtomicInteger();

//
//        for (String psmF : psmfiles) {
//         System.out.println(spectrumFolder+psmF);
//            System.out.println(spectrumFolder+psmfilesOuts[icount.get()]);
//            icount.incrementAndGet();
//
//
//        }

        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if (!Utils.bzOnlyRawdata) {
                spec.readFromeTrailFile(file, false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile, true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();


            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));


            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if (RTWidnowDistance != null && RTWidnowDistance.length > 0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
//            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName, mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection();
            iswc.IsolationWindowCollection_paralle(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName, mapRTMSoneTrail);


            //read file from top1 fdr<0.01 dataset
            for (String psmF : psmfiles) {


                ArrayList<PSMResult> arrPsmResults = SharePeaks.ReadPSMFromResultFile(spectrumFolder+psmF);


                FileWriter fileWriter = new FileWriter(spectrumFolder+psmfilesOuts[icount.get()]);
                BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);


                bufferedWriter.write("peptide\tmsonefeature\tproteinname\tdecoyTarget\tpepmass\tmsonemass\tmsonetime\tmsoneMZ" +
                        "\tmsoneZ\tBIon\tYIon\ttargetscore\tcountDP\tfdr\tsharepeaksMAX\tsharePeaksAll" +
                        "\tstrSharePeakMax\tstrSharePeakAll\tshareMs1feature\tshareMS1MZ\tshareMS1time\n");


                ProgressBar pb = new ProgressBar("Progress", arrPsmResults.size());
                pb.start();

                Map<PSMResult, MatchPeptideWithMultiCharge> mapPSMBYIon = Collections.synchronizedMap(new ConcurrentHashMap<>());//new HashMap<>();


                //the first psm
    //            arrPsmResults.get(0).sharepeaksMax = 0;
    //            MSTwoTrailSet windowTop = iswc.FindWindowWithMZ( arrPsmResults.get(0).msoneMZ, 1);
    //            double[] arrMSTwoRts = windowTop.arrMSTwoRts;
    //            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = windowTop.mapRTMStwoTrail;
    //            MatchPeptideWithMultiCharge matchHighPSM = getBYIonsMatchWithPSM(arrPsmResults.get(0), windowTop, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
    //            matchHighPSM.setMS2TrailSelected();
    //            mapPSMBYIon.put(arrPsmResults.get(0),matchHighPSM);
    //            bufferedWriter.write(arrPsmResults.get(0) +"\t0\t0\t0\n");

                arrPsmResults.stream().parallel().forEach(
                        y -> {
                            pb.step();
    //                        MSTwoTrailSet window = iswc.FindWindowWithMZ( y.msoneMZ, 1);
                            double dmatch = -100000000;
                            List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(y.msoneMZ, 1, true);
    //                        if(y.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
    //                        {
    //                            System.out.println(arrwindow.size());
    //                        }
                            for (MSTwoTrailSet window : arrwindow) {
                                double[] arrMSTwoRts = window.arrMSTwoRts;
                                Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
                                MatchPeptideWithMultiCharge matchPSM = getBYIonsMatchWithPSM(y, window, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail, arrRts);
    //                            if (matchPSM == null)
    //                                System.out.println(y.strpeptide + " " + y.msoneZ + " " + y.msonetime + " " + y.msoneMZ);
                                if (matchPSM != null) {
    //                                if(y.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
    //                                {
    //                                    System.out.println(matchPSM.dMatch);
    //                                }
                                    if (matchPSM.dMatch > dmatch) {//save the max MATCH
                                        mapPSMBYIon.put(y, matchPSM);
                                        dmatch = matchPSM.dMatch;
                                    }
                                }
                            }
                        }
                );
                pb.stop();

    //            List<PSMResult> arrPsm = Collections.synchronizedList( mapPSMBYIon.keySet().stream().sorted(Comparator.comparingDouble(PSMResult::getTargetscore).reversed()).collect(Collectors.toList()));
    //            List<PSMResult> arrPsmtest= Collections.synchronizedList(new ArrayList<>());
                for (int i = 0; i < arrPsmResults.size(); i++) {
    //            for (int i = 0; i <10000; i++) {
                    arrPsmResults.get(i).iNo = i;
    //                arrPsmtest.add(arrPsmResults.get(i));
                    PSMResult curPsm = arrPsmResults.get(i);
                    if (mapPSMBYIon.get(curPsm) != null) {
                        curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALLSpeed();
                        curPsm.strSharePeaksAllInfo = mapPSMBYIon.get(curPsm).strSharePeaksAllInfo;
                        mapPSMBYIon.get(curPsm).setMS2TrailSelectedSpeed();


                    }
                }

                if(Utils.bzOutputSharePeaksMax) {


                    ProgressBar pbafter = new ProgressBar("Progress sharePeaksMax", arrPsmResults.size());
                    pbafter.start();

                    //            for (int i = 0; i < arrPsm.size(); i++) {
                    //            arrPsm.stream().parallel().forEach(x->
                    //            arrPsmResults.stream().parallel().forEach(x->
                    arrPsmResults.stream().parallel().forEach(x ->
                    {


                        pbafter.step();
                        //                if(i%100000==0)
                        //                    System.out.println(i+"/"+arrPsm.size());
                        //get rt range[-1min,1min] and the same mz window from [0,i-1]
                        PSMResult curPsm = x;

                        int iSharePeak = 0;
                        double maxShareMS1MZ = 0.0;
                        double maxShareMS1Time = 0.0;
//                    int iCurSharePeak = 0;
                        String strshareMax = "";
                        Long ims1id = 0L;
                        mapPSMBYIon.get(curPsm).setMS2TrailSelectedSpeed();

                        for (int is = 0; is < curPsm.iNo; is++) {

                            //                    PSMResult highPsm = arrPsmResults.get(is);
                            PSMResult highPsm = arrPsmResults.get(is);
                            if (curPsm.msonetime > highPsm.msonetime - 1 && curPsm.msonetime < highPsm.msonetime + 1
                                    && iswc.isSmaeMS1Window(curPsm.msoneMZ, highPsm.msoneMZ)) {

//                            if(iCurSharePeak>0)
//                                mapPSMBYIon.get(curPsm).setMS2TrailSelectedSpeed();
                                //                        MSTwoTrailSet window = iswc.FindWindowWithMZ(curPsm.msoneMZ, 1);

                                //                        int iCurSharePeak = getSharePeak(curPsm, highPsm,window,mapRTMSoneTrail,mapPSMBYIon,arrRts);
//                            iCurSharePeak = getSharePeakWithNumber(curPsm, highPsm, mapPSMBYIon);
                                int iCurSharePeak = getSharePeak(curPsm, highPsm, mapPSMBYIon);
                                if (iCurSharePeak > iSharePeak) {
                                    iSharePeak = iCurSharePeak;
                                    maxShareMS1MZ = highPsm.msoneMZ;
                                    maxShareMS1Time = highPsm.msonetime;
                                    ims1id = highPsm.msonefeature;
                                    strshareMax = mapPSMBYIon.get(curPsm).strSharePeaksMaxInfo;
                                }
                            }
                        }
                        curPsm.sharepeaksMax = iSharePeak;
                        curPsm.strSharePeaksMaxInfo = strshareMax;
                        curPsm.maxShareMS1MZ = maxShareMS1MZ;
                        curPsm.maxShareMS1Time = maxShareMS1Time;
                        curPsm.ims1id = ims1id;


                    });
                    pbafter.stop();
                }
                ProgressBar pbafter1 = new ProgressBar("Progress output", arrPsmResults.size());
                pbafter1.start();
    //            for (int i = 0; i < arrPsmResults.size(); i++) {
                for (int i = 0; i < arrPsmResults.size(); i++) {
                    pbafter1.step();

                    //get rt range[-1min,1min] and the same mz window from [0,i-1]
    //                PSMResult curPsm = arrPsmResults.get(i);
                    PSMResult curPsm = arrPsmResults.get(i);


                    bufferedWriter.write(curPsm + "\t" + curPsm.ims1id + "\t" + curPsm.maxShareMS1MZ + "\t" + curPsm.maxShareMS1Time + "\n");


                }
                pbafter1.stop();
                    icount.incrementAndGet();
                bufferedWriter.flush();
                bufferedWriter.close();
                fileWriter.close();
            }

        }




    }

    public void match_FDRLT01_CheckSharePeaksAtRT1MinRangeMT_CorssWindow_parallel(String mgfInfile,String rawfile, String fastaInfile,String resultPSMFileName, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();



        //read file from top1 fdr<0.01 dataset
        ArrayList<PSMResult> arrPsmResults = SharePeaks.ReadPSMFromResultFile(resultPSMFileName);


        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);



        bufferedWriter.write("peptide\tmsonefeature\tproteinname\tdecoyTarget\tpepmass\tmsonemass\tmsonetime\tmsoneMZ" +
                "\tmsoneZ\tBIon\tYIon\ttargetscore\tcountDP\tfdr\tsharepeaksMAX\tsharePeaksAll" +
                "\tstrSharePeakMax\tstrSharePeakAll\tshareMs1feature\tshareMS1MZ\tshareMS1time\n");



        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if (!Utils.bzOnlyRawdata) {
                spec.readFromeTrailFile(file, false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile, true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();


            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));


            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if (RTWidnowDistance != null && RTWidnowDistance.length > 0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
//            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName, mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection();
            iswc.IsolationWindowCollection_paralle(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);


            ProgressBar pb = new ProgressBar("Progress", arrPsmResults.size());
            pb.start();

            Map<PSMResult,MatchPeptideWithMultiCharge> mapPSMBYIon = Collections.synchronizedMap( new ConcurrentHashMap<>());//new HashMap<>();


            //the first psm
//            arrPsmResults.get(0).sharepeaksMax = 0;
//            MSTwoTrailSet windowTop = iswc.FindWindowWithMZ( arrPsmResults.get(0).msoneMZ, 1);
//            double[] arrMSTwoRts = windowTop.arrMSTwoRts;
//            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = windowTop.mapRTMStwoTrail;
//            MatchPeptideWithMultiCharge matchHighPSM = getBYIonsMatchWithPSM(arrPsmResults.get(0), windowTop, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
//            matchHighPSM.setMS2TrailSelected();
//            mapPSMBYIon.put(arrPsmResults.get(0),matchHighPSM);
//            bufferedWriter.write(arrPsmResults.get(0) +"\t0\t0\t0\n");

            arrPsmResults.stream().parallel().forEach(
                    y -> {
                        pb.step();
//                        MSTwoTrailSet window = iswc.FindWindowWithMZ( y.msoneMZ, 1);
                        double dmatch = -100000000;
                        List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(y.msoneMZ, 1,true);
//                        if(y.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
//                        {
//                            System.out.println(arrwindow.size());
//                        }
                        for(MSTwoTrailSet  window:arrwindow) {
                            double[] arrMSTwoRts = window.arrMSTwoRts;
                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
                            MatchPeptideWithMultiCharge matchPSM = getBYIonsMatchWithPSM(y, window, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail, arrRts);
//                            if (matchPSM == null)
//                                System.out.println(y.strpeptide + " " + y.msoneZ + " " + y.msonetime + " " + y.msoneMZ);
                            if(matchPSM != null) {
//                                if(y.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
//                                {
//                                    System.out.println(matchPSM.dMatch);
//                                }
                                if(matchPSM.dMatch>dmatch) {//save the max MATCH
                                    mapPSMBYIon.put(y, matchPSM);
                                    dmatch = matchPSM.dMatch;
                                }
                            }
                        }
                    }
            );
            pb.stop();

//            List<PSMResult> arrPsm = Collections.synchronizedList( mapPSMBYIon.keySet().stream().sorted(Comparator.comparingDouble(PSMResult::getTargetscore).reversed()).collect(Collectors.toList()));
//            List<PSMResult> arrPsmtest= Collections.synchronizedList(new ArrayList<>());
            for (int i = 0; i < arrPsmResults.size(); i++) {
//            for (int i = 0; i <10000; i++) {
                arrPsmResults.get(i).iNo = i;
//                arrPsmtest.add(arrPsmResults.get(i));
                PSMResult curPsm = arrPsmResults.get(i);
                if(mapPSMBYIon.get(curPsm)!=null) {
                    curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALLSpeed();
                    curPsm.strSharePeaksAllInfo = mapPSMBYIon.get(curPsm).strSharePeaksAllInfo;
                    mapPSMBYIon.get(curPsm).setMS2TrailSelectedSpeed();


                }
            }
            ProgressBar pbafter = new ProgressBar("Progress sharePeaksMax", arrPsmResults.size());
            pbafter.start();

//            for (int i = 0; i < arrPsm.size(); i++) {
//            arrPsm.stream().parallel().forEach(x->
//            arrPsmResults.stream().parallel().forEach(x->
            arrPsmResults.stream().parallel().forEach(x->
            {


                pbafter.step();
//                if(i%100000==0)
//                    System.out.println(i+"/"+arrPsm.size());
                //get rt range[-1min,1min] and the same mz window from [0,i-1]
                PSMResult curPsm = x;

                int iSharePeak = 0;
                double maxShareMS1MZ = 0.0;
                double maxShareMS1Time = 0.0;
                String strshareMax="";
                Long ims1id =0L;
                for (int is = 0; is < curPsm.iNo; is++) {

//                    PSMResult highPsm = arrPsmResults.get(is);
                    PSMResult highPsm = arrPsmResults.get(is);
                    if (curPsm.msonetime > highPsm.msonetime - 1 && curPsm.msonetime < highPsm.msonetime + 1
                            && iswc.isSmaeMS1Window(curPsm.msoneMZ, highPsm.msoneMZ)) {

//                        MSTwoTrailSet window = iswc.FindWindowWithMZ(curPsm.msoneMZ, 1);

//                        int iCurSharePeak = getSharePeak(curPsm, highPsm,window,mapRTMSoneTrail,mapPSMBYIon,arrRts);
                        int iCurSharePeak = getSharePeak(curPsm, highPsm,mapPSMBYIon);
                        if (iCurSharePeak > iSharePeak)
                        {
                            iSharePeak = iCurSharePeak;
                            maxShareMS1MZ = highPsm.msoneMZ;
                            maxShareMS1Time = highPsm.msonetime;
                            ims1id = highPsm.msonefeature;
                            strshareMax = mapPSMBYIon.get(curPsm).strSharePeaksMaxInfo;
                        }
                    }


                }
                curPsm.sharepeaksMax = iSharePeak;
                curPsm.strSharePeaksMaxInfo = strshareMax;
                curPsm.maxShareMS1MZ = maxShareMS1MZ;
                curPsm.maxShareMS1Time=maxShareMS1Time;
                curPsm.ims1id=ims1id;




            });
            pbafter.stop();
            ProgressBar pbafter1 = new ProgressBar("Progress output", arrPsmResults.size());
            pbafter1.start();
//            for (int i = 0; i < arrPsmResults.size(); i++) {
            for (int i = 0; i < arrPsmResults.size(); i++) {
                pbafter1.step();

                //get rt range[-1min,1min] and the same mz window from [0,i-1]
//                PSMResult curPsm = arrPsmResults.get(i);
                PSMResult curPsm = arrPsmResults.get(i);


                bufferedWriter.write(curPsm +"\t"+curPsm.ims1id+"\t"+curPsm.maxShareMS1MZ+"\t"+curPsm.maxShareMS1Time+"\n");


            }
            pbafter1.stop();



        }


        bufferedWriter.flush();
        bufferedWriter.close();

    }

    public void match_FDRLT01_CheckSharePeaksAtRT1MinRange(String mgfInfile,String rawfile, String fastaInfile,String resultPSMFileName, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();



        //read file from top1 fdr<0.01 dataset
        ArrayList<PSMResult> arrPsmResults = SharePeaks.ReadPSMFromResultFile(resultPSMFileName);

        // read fasta
//        ArrayList<Genome> genomes = FastaFile.ReadFile(fastaInfile);
//        ArrayList<Genome> genomes = FastaFile.ReadFormPSMList(arrPsmResults);
//
//
//        // get unique set of peptides:
//        Collections.reverse(genomes); // put real ones on top
//        List<Peptide> peps = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
//        double[] pepsMass = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().mapToDouble(x -> x.mass).toArray();
//
//        System.out.println("The number of peps is: " + peps.size());




        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);



        bufferedWriter.write("peptide\tmsonefeature\tproteinname\tdecoyTarget\tpepmass\tmsonemass\tmsonetime\tmsoneMZ" +
                "\tmsoneZ\tBIon\tYIon\ttargetscore\tcountDP\tfdr\tsharepeaksMAX\tsharePeaksAll" +
                "\tstrSharePeakMax\tstrSharePeakAll\tshareMs1feature\tshareMS1MZ\tshareMS1time\n");



        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if (!Utils.bzOnlyRawdata) {
                spec.readFromeTrailFile(file, false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile, true);//读取raw的isotope文件
            spec.clearRTSandMapRT();


//            icount = spec.arrMSOneTrail.size();
            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();


            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));


            System.out.println("The number of precursor is: " + spec.arrMSOneTrail.size());


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if (RTWidnowDistance != null && RTWidnowDistance.length > 0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

            //             IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges.out", mapRTMSoneTrail);


            ProgressBar pb = new ProgressBar("Progress", arrRts.length);
            pb.start();

            Map<PSMResult,MatchPeptideWithMultiCharge> mapPSMBYIon = new HashMap<>();

            //the first psm
            arrPsmResults.get(0).sharepeaksMax = 0;
            MSTwoTrailSet windowTop = iswc.FindWindowWithMZ( arrPsmResults.get(0).msoneMZ, 1);
            double[] arrMSTwoRts = windowTop.arrMSTwoRts;
            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = windowTop.mapRTMStwoTrail;
            MatchPeptideWithMultiCharge matchHighPSM = getBYIonsMatchWithPSM(arrPsmResults.get(0), windowTop, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
            matchHighPSM.setMS2TrailSelected();
            mapPSMBYIon.put(arrPsmResults.get(0),matchHighPSM);
            bufferedWriter.write(arrPsmResults.get(0) +"\t0\t0\t0\n");

            for (int i = 1; i < arrPsmResults.size(); i++) {
                //get rt range[-1min,1min] and the same mz window from [0,i-1]
                PSMResult curPsm = arrPsmResults.get(i);

                int iSharePeak = 0;
                double maxShareMS1MZ = 0.0;
                double maxShareMS1Time = 0.0;
                String strshareMax="";
                Long ims1id =0L;
                for (int is = 0; is < i; is++) {

                    PSMResult highPsm = arrPsmResults.get(is);
                    if (curPsm.msonetime > highPsm.msonetime - 1 && curPsm.msonetime < highPsm.msonetime + 1
                            && iswc.isSmaeMS1Window(curPsm.msoneMZ, highPsm.msoneMZ)) {
                        MSTwoTrailSet window = iswc.FindWindowWithMZ(curPsm.msoneMZ, 1);

                        int iCurSharePeak = getSharePeak(curPsm, highPsm,window,mapRTMSoneTrail,mapPSMBYIon,arrRts);
                        if (iCurSharePeak > iSharePeak)
                        {
                            iSharePeak = iCurSharePeak;
                            maxShareMS1MZ = highPsm.msoneMZ;
                            maxShareMS1Time = highPsm.msonetime;
                            ims1id = highPsm.msonefeature;
                            strshareMax = mapPSMBYIon.get(curPsm).strSharePeaksMaxInfo;
                        }
                    }


                }
                curPsm.sharepeaksMax = iSharePeak;
                curPsm.strSharePeaksMaxInfo = strshareMax;
                if(mapPSMBYIon.get(curPsm)!=null) {
//                    curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALL();
                    curPsm.sharePeaksAll = mapPSMBYIon.get(curPsm).getSharePeakALLSpeed();
                    curPsm.strSharePeaksAllInfo = mapPSMBYIon.get(curPsm).strSharePeaksAllInfo;
                    mapPSMBYIon.get(curPsm).setMS2TrailSelected();


                }else {
                    MSTwoTrailSet windowCurPSM = iswc.FindWindowWithMZ( curPsm.msoneMZ, 1);
                    double[] arrMSTwoRtsCurPSM = windowCurPSM.arrMSTwoRts;
                    Map<Double, List<MSTwoTrail>> mapRTMStwoTrailCurPSM = windowCurPSM.mapRTMStwoTrail;
                    MatchPeptideWithMultiCharge matchHighPSMCurPSM = getBYIonsMatchWithPSM(curPsm, windowCurPSM, mapRTMSoneTrail, arrMSTwoRtsCurPSM, mapRTMStwoTrailCurPSM,arrRts);
                    matchHighPSMCurPSM.setMS2TrailSelected();
                    mapPSMBYIon.put(curPsm,matchHighPSMCurPSM);
                }

                bufferedWriter.write(curPsm +"\t"+ims1id+"\t"+maxShareMS1MZ+"\t"+maxShareMS1Time+"\n");


            }

            pb.stop();

        }


        bufferedWriter.flush();
        bufferedWriter.close();

    }

    private int getSharePeak(PSMResult curPsm, PSMResult highPsm, MSTwoTrailSet window, TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail, Map<PSMResult, MatchPeptideWithMultiCharge> mapPSMBYIon, double[] arrRts) {
        int iPeakSahre = 0;
        if (window != null) {
            double[] arrMSTwoRts = window.arrMSTwoRts;
            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;

            // recode each b y ion for best match ion


//            MSOneTrail msonetrail = new MSOneTrail( -1,curPsm.msoneMZ,curPsm.msonetime,curPsm.msoneZ,0,0,0.0);

            MatchPeptideWithMultiCharge matchHighPSM = mapPSMBYIon.get(highPsm);
            if(mapPSMBYIon.get(highPsm) ==null) {
                matchHighPSM= getBYIonsMatchWithPSM(highPsm, window, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
//                matchHighPSM.setMS2TrailSelected();
                mapPSMBYIon.put(highPsm,matchHighPSM);

            }

            MatchPeptideWithMultiCharge matchCurPSM = mapPSMBYIon.get(curPsm);
            if(mapPSMBYIon.get(curPsm) ==null) {
                matchCurPSM= getBYIonsMatchWithPSM(curPsm, window, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail,arrRts);
//                matchCurPSM.setMS2TrailSelected();
                mapPSMBYIon.put(curPsm,matchCurPSM);

            }
//            MatchPeptideWithMultiCharge matchCurPSM = getBYIonsMatchWithPSM(curPsm, window, mapRTMSoneTrail, arrMSTwoRts, mapRTMStwoTrail);

            if (matchCurPSM!=null) {
                iPeakSahre = matchCurPSM.getSharePeakWithOtherMatchPSMSpeed(matchHighPSM);
            }
//                iPeakSahre = matchCurPSM.getSharePeakWithOtherMatchPSM(matchHighPSM);
            else
            {
                System.out.println("msontrail not match:"+ curPsm.msonefeature+" "+curPsm.msoneMZ+" "+ curPsm.msonetime);

            }




        }//end of mz of peptide
        return iPeakSahre;
    }

    private int getSharePeakWithNumber(PSMResult curPsm, PSMResult highPsm,  Map<PSMResult, MatchPeptideWithMultiCharge> mapPSMBYIon) {
        int iPeakSahre = 0;

        MatchPeptideWithMultiCharge matchHighPSM = mapPSMBYIon.get(highPsm);
        if(mapPSMBYIon.get(highPsm) ==null) {
            System.out.println(highPsm.strpeptide+" "+highPsm.msoneMZ+" "+highPsm.msonetime+"No psm Match");

        }

        matchHighPSM.setMS2TrailMAXSelectedSpeed();

        MatchPeptideWithMultiCharge matchCurPSM = mapPSMBYIon.get(curPsm);
        if(mapPSMBYIon.get(curPsm) ==null) {
            System.out.println(curPsm.strpeptide+" "+curPsm.msoneMZ+" "+curPsm.msonetime+"No psm Match");

        }


        if (matchCurPSM!=null) {
            iPeakSahre = matchCurPSM.getSharePeakWithOtherMatchWithNumberPSMSpeed();
        }
//                iPeakSahre = matchCurPSM.getSharePeakWithOtherMatchPSM(matchHighPSM);
        else
        {
            System.out.println("msontrail not match:"+ curPsm.msonefeature+" "+curPsm.msoneMZ+" "+ curPsm.msonetime);

        }




        return iPeakSahre;
    }
    private int getSharePeak(PSMResult curPsm, PSMResult highPsm,  Map<PSMResult, MatchPeptideWithMultiCharge> mapPSMBYIon) {
        int iPeakSahre = 0;

            MatchPeptideWithMultiCharge matchHighPSM = mapPSMBYIon.get(highPsm);
            if(mapPSMBYIon.get(highPsm) ==null) {
                System.out.println(highPsm.strpeptide+" "+highPsm.msoneMZ+" "+highPsm.msonetime+"No psm Match");

            }

            MatchPeptideWithMultiCharge matchCurPSM = mapPSMBYIon.get(curPsm);
            if(mapPSMBYIon.get(curPsm) ==null) {
                System.out.println(curPsm.strpeptide+" "+curPsm.msoneMZ+" "+curPsm.msonetime+"No psm Match");

            }


            if (matchCurPSM!=null) {
                iPeakSahre = matchCurPSM.getSharePeakWithOtherMatchPSMSpeed(matchHighPSM);
            }
//                iPeakSahre = matchCurPSM.getSharePeakWithOtherMatchPSM(matchHighPSM);
            else
            {
                System.out.println("msontrail not match:"+ curPsm.msonefeature+" "+curPsm.msoneMZ+" "+ curPsm.msonetime);

            }




        return iPeakSahre;
    }
    private static MatchPeptideWithMultiCharge getBYIonsMatchWithPSM(PSMResult curPsm, MSTwoTrailSet window, TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail,
                                                                     double[] arrMSTwoRts, Map<Double, List<MSTwoTrail>> mapRTMStwoTrail, double[] arrRts) {

        double dPeptideMatchMax = -1000000.0;

        MSTwoTrail[] curBionMatchTrail = new MSTwoTrail[curPsm.strpeptide.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];
        MSTwoTrail[] curYionMatchTrail = new MSTwoTrail[curPsm.strpeptide.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];



        int[] arrbLossCount = new int[3];
        int[] arryLossCount = new int[3];
        MatchPeptideWithMultiCharge matchPeptdideMax = null;
        double diffrt = 10000.0;
        //find the ms1 rt
        int iMS1RTPosBegin = Arrays.binarySearch(arrRts, curPsm.msonetime - Utils.thresholdRT);
        iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
        if (iMS1RTPosBegin > 0) {
            iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position
        }

        //get the right pos
        int irp=iMS1RTPosBegin;
        for (;
             iMS1RTPosBegin < arrRts.length && arrRts[iMS1RTPosBegin] < curPsm.msonetime + Utils.thresholdRT;//2.0/60; is a threshold,iMS1RTPosBegin-1 is right  pos
             iMS1RTPosBegin++)
        {

            if(Math.abs(arrRts[iMS1RTPosBegin]-curPsm.msonetime)<diffrt)
            {
                diffrt = Math.abs(arrRts[iMS1RTPosBegin]-curPsm.msonetime);
                irp = iMS1RTPosBegin;

            }

        }


        int ipos =  Utils.binarySearch0(mapRTMSoneTrail.get(arrRts[irp]),0, mapRTMSoneTrail.get(arrRts[irp]).size(), curPsm.msoneMZ*(1- Utils.thresholdPPM/2));
        int iMS1MzPosBegin = Utils.getiPosBeginOfNearest(mapRTMSoneTrail.get(arrRts[irp]).size(), ipos);
        double diffmz=10000.0;
//        while (iMS1MzPosBegin > 0) {
//              iMS1MzPosBegin = iMS1MzPosBegin - 1;//move to less than position
//              if(iMS1MzPosBegin > 0 &&  mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getMz()> mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin-1).getMz()) break;//PEPTIDE的质量可能有相等的情况，这是peptide有序情况
//       }
        if (iMS1MzPosBegin > 0) {
            iMS1MzPosBegin = iMS1MzPosBegin - 1;//move to less than position
        }
        int imzp = iMS1MzPosBegin;
//        if(curPsm.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR")) {
//            System.out.println("condition diff in:"+ iMS1MzPosBegin +" "+ mapRTMSoneTrail.get(arrRts[irp]).size()+  (iMS1MzPosBegin < mapRTMSoneTrail.get(arrRts[irp]).size())
//                    +" "+mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getMz()
//                    +" "+ curPsm.msoneMZ*(1+ Utils.thresholdPPM/2.0)
//            );
//        }
        for (;
             (( iMS1MzPosBegin < mapRTMSoneTrail.get(arrRts[irp]).size()) && (mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getMz() < (curPsm.msoneMZ*(1+ Utils.thresholdPPM/2.0))));//2.0/60; is a threshold,iMS1RTPosBegin-1 is right  pos
             iMS1MzPosBegin++)
        {
            if(mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getZ()==curPsm.msoneZ) {

                if (Math.abs(mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getMz() - curPsm.msoneMZ) < diffmz) {

                    diffmz = Math.abs(mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getMz() - curPsm.msoneMZ);
                    imzp = iMS1MzPosBegin;
//                    if(curPsm.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR")) {
//                        System.out.println("msontrail diff in:"+ diffmz+" "+mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getZ()
//                                +" "+mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getMz()
//                                +" "+mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getStratRt()
//                                +" "+mapRTMSoneTrail.get(arrRts[irp]).get(iMS1MzPosBegin).getEndRt()
//                                +" "+imzp);
//                    }
                }
            }
        }

//        if(curPsm.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR")) {
//            int ib = imzp - 10;
//            if (ib < 0) ib = 0;
//            int ie = imzp + 10;
//            if (ie > mapRTMSoneTrail.get(arrRts[irp]).size()) ie = mapRTMSoneTrail.get(arrRts[irp]).size();
//            StringBuffer strinfo = new StringBuffer();
//            for (int it = ib; it < ie; it++) {
//                MSOneTrail msonetrailt = mapRTMSoneTrail.get(arrRts[irp]).get(it);
//                strinfo.append("ms1 info[-5,5] :" + (it - ib) + " " + msonetrailt.getId() + " rt:" + msonetrailt.getRt() + " mz:" + msonetrailt.getMz()+" "+msonetrailt.getStratRt()+" "+msonetrailt.getEndRt()+" Z:"+msonetrailt.getZ());
//            }
//            System.out.println("msontrail not in:"+ curPsm.msonefeature+" "+ curPsm.msoneMZ+" "+ curPsm.msoneZ+" "+ curPsm.msonetime+" RT"+arrRts[irp]+" MZ"
//                    +strinfo+" "+curPsm.strpeptide);
//        }

//        System.out.println("psm info rt:"+curPsm.msonetime+" mz:"+curPsm.msoneMZ);
//        System.out.println("psm info match:"+ curPsm.msonefeature+" "+curPsm.msoneMZ+" "+ curPsm.msonetime);

        MSOneTrail msonetrail = null;
        if((imzp>=0) &&(imzp<mapRTMSoneTrail.get(arrRts[irp]).size())) {
            msonetrail = mapRTMSoneTrail.get(arrRts[irp]).get(imzp);
//            if(curPsm.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR")) {
//
//                System.out.println("get ms1 info " + msonetrail.getId() + " rt:" + msonetrail.getRt() + " mz:" + msonetrail.getMz());
//            }
//            int ib = imzp-5;
//            if(ib<0) ib =0;
//            int ie = imzp+5;
//            if(ie>mapRTMSoneTrail.get(arrRts[irp]).size()) ie=mapRTMSoneTrail.get(arrRts[irp]).size();
//            for(int it=ib;it<ie ;it++)
//            {
//                MSOneTrail msonetrailt = mapRTMSoneTrail.get(arrRts[irp]).get(it);
//                System.out.println("ms1 info[-5,5] :"+(it-ib)+" "+ msonetrailt.getId()+" rt:" + msonetrailt.getRt() + " mz:" + msonetrailt.getMz());
//            }
        }
        else
        {
            msonetrail = new MSOneTrail( -1, curPsm.msoneMZ, curPsm.msonetime, curPsm.msoneZ,0,0,0.0);

             int ib = imzp-5;
            if(ib<0) {
                ib =0;
            }
            int ie = imzp+5;
            if(ie>mapRTMSoneTrail.get(arrRts[irp]).size()) {
                ie=mapRTMSoneTrail.get(arrRts[irp]).size();
            }
            StringBuffer strinfo=new StringBuffer();
            for(int it=ib;it<ie ;it++)
            {
                MSOneTrail msonetrailt = mapRTMSoneTrail.get(arrRts[irp]).get(it);
                strinfo.append("ms1 info[-5,5] :"+(it-ib)+" "+ msonetrailt.getId()+" rt:" + msonetrailt.getRt() + " mz:" + msonetrailt.getMz());
            }
            System.out.println("msontrail not in:"+ curPsm.msonefeature+" "+ curPsm.msoneMZ+" "+ curPsm.msonetime+" RT"+arrRts[irp]+" MZ"
                    +strinfo+" "+curPsm.strpeptide);
            return null;

        }

        int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getStratRt() - Utils.thresholdRT);
        iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
        if (iTwoRTPosBegin > 0) {
            iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
        }


        for (int iRTPos = iTwoRTPosBegin;
             iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getEndRt() + Utils.thresholdRT;
             iRTPos++) {


            Peptide pepsiPos = new Peptide(curPsm.proteinname, curPsm.strpeptide);

            pepsiPos.GenerateIons();
            double[][] arrIntensity = null;

//            if(curPsm.strpeptide.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
//            {
//                System.out.println( curPsm.strpeptide+" "+arrMSTwoRts[iRTPos]+" "+msonetrail.getId()+" "+msonetrail.getStratRt()+" "+msonetrail.getEndRt());
//            }

            MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiChargeFragmentLoss(pepsiPos, mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                    Config.scoreMode, null, 0, msonetrail, curBionMatchTrail, curYionMatchTrail, arrbLossCount, arryLossCount
                    ,arrIntensity);
            if (matchPep != null) {
                matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错

                matchPep.calculateCombineScore();
                matchPep.setDms2rt(arrMSTwoRts[iRTPos]);

//                System.out.println("matchPep:"+matchPep.dMatch +" "+ dPeptideMatchMax);

                if (matchPep.dMatch > dPeptideMatchMax) {
                    dPeptideMatchMax = matchPep.dMatch;
                    matchPeptdideMax = matchPep;
                }
            }

        }
        return matchPeptdideMax;
    }
    private static void MS1MatchPeptide_top6PredictIntensity(AtomicInteger icount, List<Peptide> peps, BufferedWriter bufferedWriter, IsolationWindowCollection iswc,
                                        Map<String, MS2Intensity> MapMS2Intensity, MSOneTrail msonetrail,
                                        LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide, Peptide pepsiPos, MSTwoTrailSet window) {
        if (window != null) {
            double[] arrMSTwoRts = window.arrMSTwoRts;
            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;


            double dPeptideMatchScoreSum = 0.0;//用于记录每个peptide在整个保留时间的得分
            double dPeptideMatchMax = -1000000.0;

            int iMS1TrailInMS2TrailCount = 0;

            MatchPeptideWithMultiCharge matchPeptdideMax = null;
            MatchPeptideWithMultiCharge matchPeptdideAll = null;

            // recode each b y ion for best match ion
            MSTwoTrail[] BionMatchTrail = new MSTwoTrail[pepsiPos.composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];
            MSTwoTrail[] YionMatchTrail = new MSTwoTrail[pepsiPos.composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];

            int[] arrbLossCount = new int[3];
            int[] arryLossCount = new int[3];

//            int[] arrbTop6IntensityPosCharge1 = new int[6];
//            int[] arryTop6IntensityPosCharge1 = new int[6];
//            int[] arrallTop6IntensityPosCharge1 = new int[6];


            pepsiPos.GenerateIons();
            double[][] arrIntensity = null;
            MS2Intensity ms2Int = MapMS2Intensity.get(pepsiPos.composition+ msonetrail.getZ());
            if (ms2Int!=null) {

                arrIntensity = ms2Int.arrdIntensity;

            }


            //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
            //0725add rt[getStratRt,getEndRt]
//                                                    int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getRt() - Utils.thresholdRT);
            int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getStratRt() - Utils.thresholdRT);
            iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
            if (iTwoRTPosBegin > 0) {
                iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
            }



            for (int iRTPos = iTwoRTPosBegin;
                 iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getEndRt() + Utils.thresholdRT;
                 iRTPos++) {
                //TEST_DEBUG
               /* if(DebugRTTime(msonetrail.getRt(),peps.get(iPos).composition))
                {
                    System.out.println(msonetrail.getRt()+"  "+peps.get(iPos).composition);
                    try {
                        FileWriter fileWriterTest = new FileWriter(psmOutfile+"html/"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass()+".html");
                        BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterTest);
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart1Title+"_"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass());
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi1Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));

                        bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi2Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart2Ms2trailData);
                        int lenth = mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size();
                        double[][] ddData= new double[lenth][3];
//                                                                bufferedWriterTest.write(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass()+"\n");
                        int iL =0;
                        for(MSTwoTrail mstwot:mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]))
                        {
                            ddData[iL][0]=mstwot.getMzApex();
                            ddData[iL][1]=mstwot.getPeakArea();
                            ddData[iL][2]=iL+1;
                            iL++;


                        }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().mapToDouble(MSTwoTrail::getMzApex).toArray())+"\n");

                        bufferedWriterTest.write(JSONWriter.valueToString(ddData)+"\n");




                        double  ponDA = 1.00727647;

                        double[][] ddBData= new double[peps.get(iPos).b_ions.length][3];


                        for(int ibs=0;ibs<peps.get(iPos).b_ions.length;ibs++)
                        {
                            ddBData[ibs][0]=peps.get(iPos).b_ions[ibs];
                            ddBData[ibs][1]=500;
                            ddBData[ibs][2]=ibs+1;


                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3BData);


                        bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");
                        double [] b_Charge = peps.get(iPos).b_ions.clone();
                        for(int ibs=0;ibs<b_Charge.length;ibs++)
                        {
                            ddBData[ibs][0]=(b_Charge[ibs]+ponDA)/2;//已经加了一个
                            ddBData[ibs][1]=400;
                            ddBData[ibs][2]=ibs+1;
//                                                                    b_Charge[ibs] = (b_Charge[ibs]+2*ponDA)/2;
                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3B2ChargeData);

                        bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");

//                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)//3charge
//                                                                {
//                                                                    ddBData[ibs][0]=(ddBData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddBData[ibs][1]=300;
//                                                                    ddBData[ibs][2]=ibs+1;
////                                                                    b_Charge[ibs] = (b_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");


                        double[][] ddYData= new double[peps.get(iPos).y_ions.length][3];


                        for(int ibs=0;ibs<peps.get(iPos).y_ions.length;ibs++)
                        {
                            ddYData[ibs][0]=peps.get(iPos).y_ions[ibs];
                            ddYData[ibs][1]=500;
//                                                                    ddYData[ibs][2]=ibs+1;
                            ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3YData);

                        bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                        double [] y_Charge = peps.get(iPos).y_ions.clone();
                        for(int ibs=0;ibs<y_Charge.length;ibs++)
                        {
                            ddYData[ibs][0]=(y_Charge[ibs]+ponDA)/2;//已经加了一个
                            ddYData[ibs][1]=400;
//                                                                    ddYData[ibs][2]=ibs+1;
                            ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


//                                                                    y_Charge[ibs]=(y_Charge[ibs]+2*ponDA)/2;
                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3Y2ChargeData);

                        bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(y_Charge)+"\n");
//                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
//                                                                {
//                                                                    ddYData[ibs][0]=(ddYData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddYData[ibs][1]=300;
//                                                                    ddYData[ibs][2]=ibs+1;
////                                                                    y_Charge[ibs]=(y_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart4End);

                        bufferedWriterTest.flush();
                        bufferedWriterTest.close();

                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }*/

                ////为获得最大值，或者是加和值，每一次重复计算
                double score = 0.0;

//                                                       MatchPeptdide matchPep = DbMatch.PepMatch(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]), Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(), msonetrail.getId());

                //    System.out.println(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass());




//                if(pepsiPos.composition.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
//                {
//                    System.out.println( pepsiPos.composition+" "+arrMSTwoRts[iRTPos]+" "+msonetrail.getId()+" "+msonetrail.getStratRt()+" "+msonetrail.getEndRt());
//                }

                MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiChargeFragmentLoss_top6PredictIntensity(pepsiPos, mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                        Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail,arrbLossCount,arryLossCount
                        ,arrIntensity,ms2Int);
               /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                        Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
                        msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
                if (matchPep != null) {
                    matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错

//                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);

                    if (ms2Int!=null) {
                        matchPep.dPrositRT = ms2Int.dPredictRT;
                    }
                    matchPep.calculateCombineScore();
                    matchPep.setDms2rt(arrMSTwoRts[iRTPos]);


                    matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);


                    if(Utils.bzCheckMS2TrailContainMS1Trail)
                    {
                        for (int iMassPos = 0; iMassPos < msonetrail.getArrMS1Mzs().length; iMassPos++) {
                            for(int jr = 0; (window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos])!=null) && jr< window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).size(); jr++)
                            {
                                if (msonetrail.getArrMS1Mzs()[iMassPos] == window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).get(jr).getMzApex()) {
                                    iMS1TrailInMS2TrailCount ++;
                                    break;
                                }
                            }

                        }
//                                                                matchPep.setStrMS1trailInMS2TrailCount(""+iMS1TrailInMS2TrailCount);


                    }
                    dPeptideMatchScoreSum += matchPep.dMatch;//将在时间内该Peptide所有的可以匹配的值全部加起来
                    if(Utils.PeptideMatchMode==Utils.PeptideMatchEnum.COMBINEMODE) {
                        if (matchPeptdideAll == null) {
                            matchPeptdideAll = matchPep;

                        } else {
                            matchPeptdideAll.combineResult(matchPep);
//                                                                matchPeptdideAll.combineResultSetRepeatCoverCount(matchPep);//update repeat and coverd
                        }
                    }
                    if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
                    {
                        if (matchPep.dMatch > dPeptideMatchMax) {
                            dPeptideMatchMax = matchPep.dMatch;
                            matchPeptdideMax = matchPep;
                        }else
                        {
                            matchPep = null;
                        }
                    }



                    //debug
//                    if (matchPep != null && matchPep.arr_bSum[0]+matchPep.arr_ySum[0] >10 && ms2Int!=null
//                            && ms2Int.arrallTop6IntensityPos!=null){
//                        if(Utils.iDebugoutput > 0)
//                        {
//                            try {
//                                Utils.bufferedWriterDebugInfo.write("arrallTop6IntensityCharge1:"+ JSONWriter.valueToString(ms2Int.arrallTop6IntensityCharge1)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("arrallTop6IntensityCharge1Pos:"+JSONWriter.valueToString(ms2Int.arrallTop6IntensityCharge1Pos)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("arrallTop6IntensityPos:"+JSONWriter.valueToString(ms2Int.arrallTop6IntensityPos)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("arrallTop6IntensityCount:"+JSONWriter.valueToString(ms2Int.arrallTop6IntensityCount)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("arrbTop6IntensityPosCharge1:"+JSONWriter.valueToString(ms2Int.arrbTop6IntensityPosCharge1)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("arryTop6IntensityPosCharge1:"+JSONWriter.valueToString(ms2Int.arryTop6IntensityPosCharge1)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("arrdIntensity:"+JSONWriter.valueToString(ms2Int.arrdIntensity)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("bIonMatchTrail:"+JSONWriter.valueToString(BionMatchTrail)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("yIonMatchTrail:"+JSONWriter.valueToString(YionMatchTrail)+"\n");
//                                Utils.bufferedWriterDebugInfo.write("dCosWithPredictMSMSMatchedTop6AllPredict:"+matchPep+"\n");
//
//
//                            } catch (IOException e) {
//                                throw new RuntimeException(e);
//                            }
//
//
//                        }
//
//                    }


                }

////节省时间不能保存值 begin
//                                                    double score = 0.0;
//                                                    if (mapTriplePeptideRTScore[iPos]==null)
//                                                    {
////                                                        System.out.println("----"+arrMSTwoRts[iRTPos]);
//
//                                                        score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
//                                                        mapTriplePeptideRTScore[iPos] = new BitSet();
//                                                        mapTriplePeptideRTScore[iPos].set(iRTPos);
//                                                    }else
//                                                    {
//                                                        if(mapTriplePeptideRTScore[iPos].get(iRTPos)==false)
//                                                        {
//                                                            score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
////                                                            System.out.println("----~~~~~~~~~"+arrMSTwoRts[iRTPos]);
//
//                                                            mapTriplePeptideRTScore[iPos].set(iRTPos);
//
//                                                        }
//                                                    }
////节省时间不能保存值 end

//                                                    if(score>=0.4) {
//                                                        try {
//                                                            bufferedWriter.write(icount.incrementAndGet() + "\t");
//                                                            bufferedWriter.write(peps.get(iPos).id+ "\t" +peps.get(iPos).composition + "\t"+ peps.get(iPos).mass);
//                                                            bufferedWriter.write( "\t"+msonetrail.getMass());
//                                                            bufferedWriter.write("\t"+score);
//
//                                                            bufferedWriter.write("\t"+arrMSTwoRts[iRTPos]+"\n");
//                                                        } catch (IOException e) {
//                                                            e.printStackTrace();
//                                                        }
//
////                                                        System.out.println(icount.incrementAndGet() +"--------");
////                                                        System.out.println(peps.get(iPos).composition + peps.get(iPos).mass);
////                                                        System.out.println(msonetrail.getMass());
////                                                        System.out.println(score);
//
//                                                        System.out.println(arrMSTwoRts[iRTPos]);
//                                                    }
            }
            MatchPeptideWithMultiCharge mp = null;
//                                                    MatchPeptdide mp = null;
            double dPm = 0.0;
            if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
            {
                mp = matchPeptdideMax;
                dPm = dPeptideMatchMax;
                if(mp!=null)
                {
//                                                            mp.setRepeatCoverCount(matchPeptdideAll);
                    if(Utils.bzCheckMS2TrailContainMS1Trail) {
                        mp.setStrMS1trailInMS2TrailCount("" + iMS1TrailInMS2TrailCount);
                    }
                }/*else
                {
                    //TEST BUG
                    System.out.println(matchPeptdideAll);

                }*/

            } else if(Utils.PeptideMatchMode==Utils.PeptideMatchEnum.COMBINEMODE) {
                mp = matchPeptdideAll;
                dPm = dPeptideMatchScoreSum;
            }


            //peptide search mode: the more ion covered, the more score have
            if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.CoverMax) {

                if (mp != null) {
                    int length = mp.lBIonCovered.stream().distinct().toArray().length +
                            mp.lYIonCovered.stream().distinct().toArray().length;
//                                                        mp.dMatch = mp.dMatch * length;
                    mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                }

            } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.Repeat) {
                if (mp != null) {
                    int length = mp.iBRepeatCount + mp.iYRepeatCount;
                    mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                }

            } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.RepeatPlusCoverMax) {
                if (mp != null) {
                    int length = mp.iBRepeatCount + mp.iYRepeatCount + mp.lBIonCovered.stream().distinct().toArray().length +
                            mp.lYIonCovered.stream().distinct().toArray().length;

                    mp.dMatch = mp.dMatch * length;//    / mp.strPepcomposition.length();
                }

            }


            //remain candidate list
            int iCandidateListSize = Utils.iCandidateListSize;//Candidate list no
//            if (matchPeptdideAll != null) {
            if (mp != null) {
                if (listMaxiumPeptide.size() == 0) {
                    listMaxiumPeptide.add(mp);
                } else {
                    int iL = 0;
                    Boolean bzAdd = false;
                    for (; iL < iCandidateListSize && iL < listMaxiumPeptide.size(); iL++) {
                        if (listMaxiumPeptide.get(iL).dMatch < dPm) {
                            listMaxiumPeptide.add(iL, mp);//strategy one

                            bzAdd = true;
                            break;

                        }
                    }
//                                                            if (!bzAdd && iL < iCandidateListSize - 1) {
                    if (!bzAdd && iL < iCandidateListSize ) {

                        listMaxiumPeptide.add(mp);

                    }
                    if (listMaxiumPeptide.size() > iCandidateListSize) {
                        listMaxiumPeptide.removeLast();
                    }


                }
            }

        }//end of mz of peptide
    }
    private static void MS1MatchPeptide(AtomicInteger icount, List<Peptide> peps, BufferedWriter bufferedWriter, IsolationWindowCollection iswc,
                                        Map<String, MS2Intensity> MapMS2Intensity, MSOneTrail msonetrail,
                                        LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide, Peptide pepsiPos, MSTwoTrailSet window) {
        if (window != null) {
            double[] arrMSTwoRts = window.arrMSTwoRts;
            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;


            double dPeptideMatchScoreSum = 0.0;//用于记录每个peptide在整个保留时间的得分
            double dPeptideMatchMax = -1000000.0;

            int iMS1TrailInMS2TrailCount = 0;

            MatchPeptideWithMultiCharge matchPeptdideMax = null;
            MatchPeptideWithMultiCharge matchPeptdideAll = null;

            // recode each b y ion for best match ion
            MSTwoTrail[] BionMatchTrail = new MSTwoTrail[pepsiPos.composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];
            MSTwoTrail[] YionMatchTrail = new MSTwoTrail[pepsiPos.composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];
//            MSTwoTrail[] BionMatchTrail = new MSTwoTrail[peps.get(iPos).composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];
//            MSTwoTrail[] YionMatchTrail = new MSTwoTrail[peps.get(iPos).composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];

            int[] arrbLossCount = new int[3];
            int[] arryLossCount = new int[3];




            //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
            //0725add rt[getStratRt,getEndRt]
//                                                    int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getRt() - Utils.thresholdRT);
            int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getStratRt() - Utils.thresholdRT);
            iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
            if (iTwoRTPosBegin > 0) {
                iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
            }


//
//                                                    for (int iRTPos = iTwoRTPosBegin;
//                                                         iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getRt() + Utils.thresholdRT;
//                                                         iRTPos++) {

            for (int iRTPos = iTwoRTPosBegin;
                 iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getEndRt() + Utils.thresholdRT;
                 iRTPos++) {
                //TEST_DEBUG
               /* if(DebugRTTime(msonetrail.getRt(),peps.get(iPos).composition))
                {
                    System.out.println(msonetrail.getRt()+"  "+peps.get(iPos).composition);
                    try {
                        FileWriter fileWriterTest = new FileWriter(psmOutfile+"html/"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass()+".html");
                        BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterTest);
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart1Title+"_"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass());
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi1Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));

                        bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi2Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart2Ms2trailData);
                        int lenth = mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size();
                        double[][] ddData= new double[lenth][3];
//                                                                bufferedWriterTest.write(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass()+"\n");
                        int iL =0;
                        for(MSTwoTrail mstwot:mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]))
                        {
                            ddData[iL][0]=mstwot.getMzApex();
                            ddData[iL][1]=mstwot.getPeakArea();
                            ddData[iL][2]=iL+1;
                            iL++;


                        }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().mapToDouble(MSTwoTrail::getMzApex).toArray())+"\n");

                        bufferedWriterTest.write(JSONWriter.valueToString(ddData)+"\n");




                        double  ponDA = 1.00727647;

                        double[][] ddBData= new double[peps.get(iPos).b_ions.length][3];


                        for(int ibs=0;ibs<peps.get(iPos).b_ions.length;ibs++)
                        {
                            ddBData[ibs][0]=peps.get(iPos).b_ions[ibs];
                            ddBData[ibs][1]=500;
                            ddBData[ibs][2]=ibs+1;


                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3BData);


                        bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");
                        double [] b_Charge = peps.get(iPos).b_ions.clone();
                        for(int ibs=0;ibs<b_Charge.length;ibs++)
                        {
                            ddBData[ibs][0]=(b_Charge[ibs]+ponDA)/2;//已经加了一个
                            ddBData[ibs][1]=400;
                            ddBData[ibs][2]=ibs+1;
//                                                                    b_Charge[ibs] = (b_Charge[ibs]+2*ponDA)/2;
                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3B2ChargeData);

                        bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");

//                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)//3charge
//                                                                {
//                                                                    ddBData[ibs][0]=(ddBData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddBData[ibs][1]=300;
//                                                                    ddBData[ibs][2]=ibs+1;
////                                                                    b_Charge[ibs] = (b_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");


                        double[][] ddYData= new double[peps.get(iPos).y_ions.length][3];


                        for(int ibs=0;ibs<peps.get(iPos).y_ions.length;ibs++)
                        {
                            ddYData[ibs][0]=peps.get(iPos).y_ions[ibs];
                            ddYData[ibs][1]=500;
//                                                                    ddYData[ibs][2]=ibs+1;
                            ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3YData);

                        bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                        double [] y_Charge = peps.get(iPos).y_ions.clone();
                        for(int ibs=0;ibs<y_Charge.length;ibs++)
                        {
                            ddYData[ibs][0]=(y_Charge[ibs]+ponDA)/2;//已经加了一个
                            ddYData[ibs][1]=400;
//                                                                    ddYData[ibs][2]=ibs+1;
                            ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


//                                                                    y_Charge[ibs]=(y_Charge[ibs]+2*ponDA)/2;
                        }
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart3Y2ChargeData);

                        bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(y_Charge)+"\n");
//                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
//                                                                {
//                                                                    ddYData[ibs][0]=(ddYData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddYData[ibs][1]=300;
//                                                                    ddYData[ibs][2]=ibs+1;
////                                                                    y_Charge[ibs]=(y_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                        bufferedWriterTest.write(Utils.strHtmlTemplatePart4End);

                        bufferedWriterTest.flush();
                        bufferedWriterTest.close();

                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }*/

                ////为获得最大值，或者是加和值，每一次重复计算
                double score = 0.0;

//                                                       MatchPeptdide matchPep = DbMatch.PepMatch(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]), Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(), msonetrail.getId());

                //    System.out.println(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass());

                pepsiPos.GenerateIons();
                double[][] arrIntensity = null;
                MS2Intensity ms2Int = MapMS2Intensity.get(pepsiPos.composition+ msonetrail.getZ());
                if (ms2Int!=null) {
                    arrIntensity = ms2Int.arrdIntensity;
                }


//                if(pepsiPos.composition.equals("TKDNITIYKQALLTGFNTVSVQVEMSWANLDR"))
//                {
//                    System.out.println( pepsiPos.composition+" "+arrMSTwoRts[iRTPos]+" "+msonetrail.getId()+" "+msonetrail.getStratRt()+" "+msonetrail.getEndRt());
//                }

                MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiChargeFragmentLoss(pepsiPos, mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                        Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail,arrbLossCount,arryLossCount
                        ,arrIntensity);
               /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                        Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
                        msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
                if (matchPep != null) {
                    matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错

//                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);

                    if (ms2Int!=null) {
                        matchPep.dPrositRT = ms2Int.dPredictRT;
                    }
                    matchPep.calculateCombineScore();
                    matchPep.setDms2rt(arrMSTwoRts[iRTPos]);


                    matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);


                    if(Utils.bzCheckMS2TrailContainMS1Trail)
                    {
                        for (int iMassPos = 0; iMassPos < msonetrail.getArrMS1Mzs().length; iMassPos++) {
                            for(int jr = 0; (window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos])!=null) && jr< window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).size(); jr++)
                            {
                                if (msonetrail.getArrMS1Mzs()[iMassPos] == window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).get(jr).getMzApex()) {
                                    iMS1TrailInMS2TrailCount ++;
                                    break;
                                }
                            }

                        }
//                                                                matchPep.setStrMS1trailInMS2TrailCount(""+iMS1TrailInMS2TrailCount);


                    }
                    dPeptideMatchScoreSum += matchPep.dMatch;//将在时间内该Peptide所有的可以匹配的值全部加起来
                    if (matchPeptdideAll == null) {
                        matchPeptdideAll = matchPep;

                    } else {
                        matchPeptdideAll.combineResult(matchPep);
//                                                                matchPeptdideAll.combineResultSetRepeatCoverCount(matchPep);//update repeat and coverd
                    }
                    if (matchPep.dMatch > dPeptideMatchMax) {
                        dPeptideMatchMax = matchPep.dMatch;
                        matchPeptdideMax = matchPep;
                    }
                }

////节省时间不能保存值 begin
//                                                    double score = 0.0;
//                                                    if (mapTriplePeptideRTScore[iPos]==null)
//                                                    {
////                                                        System.out.println("----"+arrMSTwoRts[iRTPos]);
//
//                                                        score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
//                                                        mapTriplePeptideRTScore[iPos] = new BitSet();
//                                                        mapTriplePeptideRTScore[iPos].set(iRTPos);
//                                                    }else
//                                                    {
//                                                        if(mapTriplePeptideRTScore[iPos].get(iRTPos)==false)
//                                                        {
//                                                            score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
////                                                            System.out.println("----~~~~~~~~~"+arrMSTwoRts[iRTPos]);
//
//                                                            mapTriplePeptideRTScore[iPos].set(iRTPos);
//
//                                                        }
//                                                    }
////节省时间不能保存值 end

//                                                    if(score>=0.4) {
//                                                        try {
//                                                            bufferedWriter.write(icount.incrementAndGet() + "\t");
//                                                            bufferedWriter.write(peps.get(iPos).id+ "\t" +peps.get(iPos).composition + "\t"+ peps.get(iPos).mass);
//                                                            bufferedWriter.write( "\t"+msonetrail.getMass());
//                                                            bufferedWriter.write("\t"+score);
//
//                                                            bufferedWriter.write("\t"+arrMSTwoRts[iRTPos]+"\n");
//                                                        } catch (IOException e) {
//                                                            e.printStackTrace();
//                                                        }
//
////                                                        System.out.println(icount.incrementAndGet() +"--------");
////                                                        System.out.println(peps.get(iPos).composition + peps.get(iPos).mass);
////                                                        System.out.println(msonetrail.getMass());
////                                                        System.out.println(score);
//
//                                                        System.out.println(arrMSTwoRts[iRTPos]);
//                                                    }
            }
            MatchPeptideWithMultiCharge mp = null;
//                                                    MatchPeptdide mp = null;
            double dPm;
            if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
            {
                mp = matchPeptdideMax;
                dPm = dPeptideMatchMax;
                if(mp!=null)
                {
//                                                            mp.setRepeatCoverCount(matchPeptdideAll);
                    if(Utils.bzCheckMS2TrailContainMS1Trail) {
                        mp.setStrMS1trailInMS2TrailCount("" + iMS1TrailInMS2TrailCount);
                    }
                }/*else
                {
                    //TEST BUG
                    System.out.println(matchPeptdideAll);

                }*/

            } else {
                mp = matchPeptdideAll;
                dPm = dPeptideMatchScoreSum;
            }


            //peptide search mode: the more ion covered, the more score have
            if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.CoverMax) {

                if (mp != null) {
                    int length = mp.lBIonCovered.stream().distinct().toArray().length +
                            mp.lYIonCovered.stream().distinct().toArray().length;
//                                                        mp.dMatch = mp.dMatch * length;
                    mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                }

            } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.Repeat) {
                if (mp != null) {
                    int length = mp.iBRepeatCount + mp.iYRepeatCount;
                    mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                }

            } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.RepeatPlusCoverMax) {
                if (mp != null) {
                    int length = mp.iBRepeatCount + mp.iYRepeatCount + mp.lBIonCovered.stream().distinct().toArray().length +
                            mp.lYIonCovered.stream().distinct().toArray().length;

                    mp.dMatch = mp.dMatch * length;//    / mp.strPepcomposition.length();
                }

            }


            //remain candidate list
            int iCandidateListSize = Utils.iCandidateListSize;//Candidate list no
            if (matchPeptdideAll != null) {
                if (listMaxiumPeptide.size() == 0) {
                    listMaxiumPeptide.add(mp);
                } else {
                    int iL = 0;
                    Boolean bzAdd = false;
                    for (; iL < iCandidateListSize && iL < listMaxiumPeptide.size(); iL++) {
                        if (listMaxiumPeptide.get(iL).dMatch < dPm) {
                            listMaxiumPeptide.add(iL, mp);//strategy one

                            bzAdd = true;
                            break;

                        }
                    }
//                                                            if (!bzAdd && iL < iCandidateListSize - 1) {
                    if (!bzAdd && iL < iCandidateListSize ) {

                        listMaxiumPeptide.add(mp);

                    }
                    if (listMaxiumPeptide.size() > iCandidateListSize) {
                        listMaxiumPeptide.removeLast();
                    }


                }
            }

        }//end of mz of peptide
    }

    public void match_TrailWithMultiRT_FitRT(String mgfInfile,String rawfile, String fastaInfile, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();

        // read fasta
        ArrayList<Genome> genomes = FastaFile.ReadFile(fastaInfile);
//        ArrayList<Genome> genomes = FastaFile.ReadFileFormPeptide(fastaInfile);


        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top
        List<Peptide> peps = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
        double[] pepsMass = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());

//        FileWriter fileWriterPeptideInfo = new FileWriter("peptideInfo.csv");
//        BufferedWriter bufferedWriterPeptideInfo = new BufferedWriter(fileWriterPeptideInfo);
//        bufferedWriterPeptideInfo.write( "peptide\tmass\tdMutationRate"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
//                + "\n");
//        for (Peptide pep:peps)
//        {
//            bufferedWriterPeptideInfo.write(pep.composition+'\t'+pep.mass+'\t'+pep.dMutationRate+'\n');
//        }
//
//        bufferedWriterPeptideInfo.flush();
//        bufferedWriterPeptideInfo.close();


        ArrayList<DbMatch> resLst = new ArrayList<>();




//        Set<Peptide>
        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tmsoneMZ\tmsoneZ\tmsonePeakAreaLocalRank" +
                "\tms2rt\tpeptideLength\tMS1PeptideMassError\tdPeptideMutationRate\tibyComplemetaryCount\tdCosWithPredictMSMS\tdMS1WithMS2TrailCoscinSimilarity\tMS2TrailCoscinSimilarity\tBIon\tYIon");
        for(int iw =0;iw<Utils.iMS2Charge;iw++) {
            bufferedWriter.write("\tdlogBionPeakAreaSum_C"+(iw+1)+"\tdlogYionPeakAreaSum_C"+(iw+1)+"\tdBionMassErrorSum_C"+(iw+1)+"\tdYionMassErrorSum_C"+(iw+1)+"\tiBionPeakHalfSum_C"+(iw+1)+"\tiYionPeakHalfSum_C"+(iw+1)+"\tbionMatchCount_C"+(iw+1)+"\tyionMatchCount" +
                    "_C"+(iw+1)+"\tdBionRetentionTimeErrorSum_C"+(iw+1)+"\tdYionRetentionTimeErrorSum_C"+(iw+1)+"\tiBConsective_C"+(iw+1)+"\tiYConsective_C"+(iw+1)+"\tdBionCosWithPredictMSMS_C"+(iw+1)+"\tdYionCosWithPredictMSMS_C"+(iw+1));
        }
        bufferedWriter.write( "\tarrMatchedBion\tarrMatchedYion\tadjustBYIon\tSoredIndex\twindowSize\tdBionCos60WithPredictMSMS_C1\tdCosWithPredictMSMSMatchedPredict\tprositRT" +
                "\tibAllMatchedWithloss\tibMatchNormalAndH2Oloss\tibMatchNormalAndNH3loss\tiyAllMatchedWithloss\tiyMatchNormalAndH2Oloss\tiyMatchNormalAndNH3loss\tisRawdata\tcount"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
                + "\n");



//        Map<Integer,Map<Integer,Double>> mapTriplePeptideRTScore = new TreeMap<>();

//        BitSet mapTriplePeptideRTScore = new BitSet(peps.size());
//        BitSet[] mapTriplePeptideRTScore = new BitSet[peps.size()];
//        Map<Pair<Integer,Integer>,Double> mapTriplePeptideRTScore = new HashMap<>();
//        for (MSTwoTrail msTwoTrail : MSTWOts.arrMSTwoTrail) {
//            double[] rts = msTwoTrail.getRts();
//
//            for(double drt:rts) {
//                if (uniqueValues.add(drt))
//                {
//                    List<MSTwoTrail> listMSTwoTrail = new ArrayList<>();
//                    listMSTwoTrail.add(msTwoTrail);
//
//                    mapRTMSTwoTrail.put(drt,listMSTwoTrail);
//                }
//                else
//                {
//                    mapRTMSTwoTrail.get(drt).add(msTwoTrail);
//                }
//            }
//        }

        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            if(!Utils.bzOnlyRawdata)
            {
                spec.readFromeTrailFile(file,false);//读取xiangyuan的MS1 FEATURE文件
                spec.generateRTSandMapRT();

            }
            spec.readFromeTrailFile(rawfile,true);//读取raw的isotope文件
            spec.clearRTSandMapRT();



            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();

//            Map<Double, List<MSOneTrail>> mapRTMSoneTrail =  spec.arrMSOneTrail.stream().
//                    collect(Collectors.groupingBy(x-> x.getRt()));

            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));





            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();
            if(RTWidnowDistance!=null && RTWidnowDistance.length>0) {
                Utils.MS1RTSpan = RTWidnowDistance[0];
            }


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

//            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
                  IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges.out",mapRTMSoneTrail);


//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2Intensity();
            Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit();

            //GET The autort information
//            Map<String, Double> MapPeptidePredictAutoRT = Main.getAutoRTInfo();
//            List<String> strpepsAutoRT = MapPeptidePredictAutoRT.keySet().stream().collect(Collectors.toList());
//            double[] dpepsAutoRT = MapPeptidePredictAutoRT.entrySet().stream().mapToDouble(x -> x.getValue()).toArray();
//
//            MapPeptidePredictAutoRT = null; //release



            //set ms1feature peakarea rank
            int ims1 = 1;
            for(MSOneTrail ms1Trail:spec.arrMSOneTrail)
            {
                long iPeakAreaRank = 0;
                MSTwoTrailSet window = iswc.FindWindowWithMZ(ms1Trail.getMz(), 1);
                //get ms1feature RT distance between Math.abs(y-v) <= Utils.thresholdRT
                int iMS1RTPosBegin = Arrays.binarySearch(arrRts, ms1Trail.getRt() - Utils.thresholdMS1RT);
                iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
                if (iMS1RTPosBegin > 0) {
                    iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position
                }

                for (int iRTPos = iMS1RTPosBegin;
                     iRTPos < arrRts.length && arrRts[iRTPos] < ms1Trail.getRt() + Utils.thresholdMS1RT;
                     iRTPos++) {
                    if(window!=null) {
                        iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
                                x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh && x.getQuantification_peaks_area() >= ms1Trail.getQuantification_peaks_area()
                        ).count();
                    }
                }
                if(iPeakAreaRank>0)
                {
                    ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));
                }else
                {
                    ms1Trail.setdPeakAreaLocalRankInTimeSpan( 0 );
                }


/*                if(window!=null) {
                    iPeakAreaRank = spec.arrMSOneTrail.stream().filter((x) ->
                            x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh && x.getQuantification_peaks_area() >= ms1Trail.getQuantification_peaks_area()
                                    && x.getRt() >= ms1Trail.getRt() - Utils.thresholdMS1RT && x.getRt() <= ms1Trail.getRt() + Utils.thresholdMS1RT
                    ).count();
                }*/

//                System.out.println("ims1:"+ ims1++  +"----"+spec.arrMSOneTrail.size());

            }

//        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(
//                "isolationWindowRanges.txt"));
//        IsolationWindowCollection iswc = (IsolationWindowCollection) ois.readObject();

//        List<double[]> larrMSTwoRts = new ArrayList<>(iswTest.windows.size());
//        List<Map<Double, List<MSTwoTrail>>> lmapRTMStwoTrail = new ArrayList<>(iswTest.windows.size());

            // TODO: READ IN SPECTRUM RANKING FILE (Added by Caroline) Implement Rank Here


            // for each isolation window collection, match once
            //Step 1: get a msonetrail
            //Step 2: sort with RT
            //Step 3:for each msonetrail get the window
            //      check all peptide with same mass  (mz +- 30 ppm)
            //         for each mstwoTrail in window (mz in window && RT +- 5s )
            //             compare mstwotrail b y ions with peptide b y ions (mz +- 30 ppm)
            //         add peptide to msonetrail as candidate peptide with confident score
//        MSTwoTrailSet MSTWOts = new MSTwoTrailSet();
//        MSTWOts.readFromeTrailFile(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window0_ms2_trails.tsv");

//        Map<Double,List<MSTwoTrail> >mapRTMSTwoTrail = new HashMap<Double, List<MSTwoTrail>>();
            TreeSet<Double> uniqueValues = new TreeSet<>();


            //////////////////////////////////
            //for a window,initial rts and rts corresponding trails
//        double[] arrMSTwoRts = MSTWOts.arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
//                .distinct()
//                .sorted()
//                .toArray();
////        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().
////                collect(Collectors.groupingBy(x-> x.getRtApex()));
//
//
//        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
            //////////////////////////////////////

//            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = new HashMap<>();
            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = Collections.synchronizedMap( new HashMap<>());

//            Map<Long, LinkedList<MatchPeptdide>> mapMSOneFeaturePeptide = new HashMap<>();

//            AtomicInteger iPrecursorNO= new AtomicInteger(1000000);
//            AtomicInteger i = new AtomicInteger();
            ProgressBar pb = new ProgressBar("Progress", arrRts.length);
            pb.start();
//            List<Double> dRTSynlist = Collections.synchronizedList(new ArrayList<>());
//            for(double dRT:arrRts)
//                dRTSynlist.add(dRT);

//            dRTSynlist.stream().parallel().forEach(
//            Arrays.stream(arrRts).forEach(
              Arrays.stream(arrRts).parallel().forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();
                        pb.step();
                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);



                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {



//                                        if(msonetrail.bAddPrecursor)
//                                        {


//                                        if (msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass >= pepsMass[0]) //filter m1 trail which less than the least peptide mass

                                        //begin of window with mz of m1 feature
////                                        MSTwoTrailSet window = iswc.FindWindow(msonetrail.getMass(), 2);
//                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                        // if mass matches, then pep match
//                                        if (window != null) {
//                                            double[] arrMSTwoRts = window.arrMSTwoRts;
//                                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;

                                        if (msonetrail.getMass() - Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {



                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window

                                            LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide = new LinkedList<MatchPeptideWithMultiCharge>();

//                                            LinkedList<MatchPeptdide> listMaxiumPeptide = new LinkedList<MatchPeptdide>();


                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况
//                                            while (iPosBegin > 0) {
//                                                iPosBegin = iPosBegin - 1;//move to less than position
//                                                if(iPosBegin > 0 && pepsMass[iPosBegin]>pepsMass[iPosBegin-1]) break;//PEPTIDE的质量可能有相等的情况，这是peptide有序情况
//                                            }

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {



                                                //begin of window with mz of ms1 feature
                                                MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
                                                // if mass matches, then pep match
                                                if (window != null) {
                                                    double[] arrMSTwoRts = window.arrMSTwoRts;
                                                    Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;


                                                    double dPeptideMatchScoreSum = 0.0;//用于记录每个peptide在整个保留时间的得分
                                                    double dPeptideMatchMax = -1000000.0;

                                                    int iMS1TrailInMS2TrailCount = 0;

                                                    MatchPeptideWithMultiCharge matchPeptdideMax = null;
                                                    MatchPeptideWithMultiCharge matchPeptdideAll = null;

                                                    // recode each b y ion for best match ion
                                                    MSTwoTrail[] BionMatchTrail = new MSTwoTrail[peps.get(iPos).composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];
                                                    MSTwoTrail[] YionMatchTrail = new MSTwoTrail[peps.get(iPos).composition.length() * Utils.iMS2Charge * Utils.iMS2FragmentTypeCount];

                                                    int[] arrbLossCount = new int[3];
                                                    int[] arryLossCount = new int[3];




                                                    //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
                                                    //0725add rt[getStratRt,getEndRt]
//                                                    int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getRt() - Utils.thresholdRT);
                                                    int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getStratRt() - Utils.thresholdRT);
                                                    iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
                                                    if (iTwoRTPosBegin > 0) {
                                                        iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
                                                    }


//
//                                                    for (int iRTPos = iTwoRTPosBegin;
//                                                         iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getRt() + Utils.thresholdRT;
//                                                         iRTPos++) {

                                                    for (int iRTPos = iTwoRTPosBegin;
                                                         iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getEndRt() + Utils.thresholdRT;
                                                         iRTPos++) {
                                                        //TEST_DEBUG
                                                       /* if(DebugRTTime(msonetrail.getRt(),peps.get(iPos).composition))
                                                        {
                                                            System.out.println(msonetrail.getRt()+"  "+peps.get(iPos).composition);
                                                            try {
                                                                FileWriter fileWriterTest = new FileWriter(psmOutfile+"html/"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass()+".html");
                                                                BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterTest);
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Title+"_"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass());
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi1Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));

                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi2Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart2Ms2trailData);
                                                                int lenth = mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size();
                                                                double[][] ddData= new double[lenth][3];
//                                                                bufferedWriterTest.write(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass()+"\n");
                                                                int iL =0;
                                                                for(MSTwoTrail mstwot:mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]))
                                                                {
                                                                    ddData[iL][0]=mstwot.getMzApex();
                                                                    ddData[iL][1]=mstwot.getPeakArea();
                                                                    ddData[iL][2]=iL+1;
                                                                    iL++;


                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().mapToDouble(MSTwoTrail::getMzApex).toArray())+"\n");

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddData)+"\n");




                                                                double  ponDA = 1.00727647;

                                                                double[][] ddBData= new double[peps.get(iPos).b_ions.length][3];


                                                                for(int ibs=0;ibs<peps.get(iPos).b_ions.length;ibs++)
                                                                {
                                                                    ddBData[ibs][0]=peps.get(iPos).b_ions[ibs];
                                                                    ddBData[ibs][1]=500;
                                                                    ddBData[ibs][2]=ibs+1;


                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3BData);


                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");
                                                                double [] b_Charge = peps.get(iPos).b_ions.clone();
                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)
                                                                {
                                                                    ddBData[ibs][0]=(b_Charge[ibs]+ponDA)/2;//已经加了一个
                                                                    ddBData[ibs][1]=400;
                                                                    ddBData[ibs][2]=ibs+1;
//                                                                    b_Charge[ibs] = (b_Charge[ibs]+2*ponDA)/2;
                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3B2ChargeData);

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");

//                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)//3charge
//                                                                {
//                                                                    ddBData[ibs][0]=(ddBData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddBData[ibs][1]=300;
//                                                                    ddBData[ibs][2]=ibs+1;
////                                                                    b_Charge[ibs] = (b_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");


                                                                double[][] ddYData= new double[peps.get(iPos).y_ions.length][3];


                                                                for(int ibs=0;ibs<peps.get(iPos).y_ions.length;ibs++)
                                                                {
                                                                    ddYData[ibs][0]=peps.get(iPos).y_ions[ibs];
                                                                    ddYData[ibs][1]=500;
//                                                                    ddYData[ibs][2]=ibs+1;
                                                                    ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3YData);

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                                                                double [] y_Charge = peps.get(iPos).y_ions.clone();
                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
                                                                {
                                                                    ddYData[ibs][0]=(y_Charge[ibs]+ponDA)/2;//已经加了一个
                                                                    ddYData[ibs][1]=400;
//                                                                    ddYData[ibs][2]=ibs+1;
                                                                    ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


//                                                                    y_Charge[ibs]=(y_Charge[ibs]+2*ponDA)/2;
                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3Y2ChargeData);

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(y_Charge)+"\n");
//                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
//                                                                {
//                                                                    ddYData[ibs][0]=(ddYData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddYData[ibs][1]=300;
//                                                                    ddYData[ibs][2]=ibs+1;
////                                                                    y_Charge[ibs]=(y_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart4End);

                                                                bufferedWriterTest.flush();
                                                                bufferedWriterTest.close();

                                                            } catch (IOException e) {
                                                                e.printStackTrace();
                                                            }
                                                        }*/

                                                        ////为获得最大值，或者是加和值，每一次重复计算
                                                        double score = 0.0;

//                                                       MatchPeptdide matchPep = DbMatch.PepMatch(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]), Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(), msonetrail.getId());

                                                        //    System.out.println(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass());

                                                        peps.get(iPos).GenerateIons();
                                                        double[][] arrIntensity = null;
                                                        MS2Intensity ms2Int = MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ());
                                                        if (ms2Int!=null) {
                                                            arrIntensity = ms2Int.arrdIntensity;
                                                        }

                                                        MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiChargeFragmentLoss(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail,arrbLossCount,arryLossCount
                                                                ,arrIntensity);
                                                       /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
                                                                msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
                                                        if (matchPep != null) {
                                                            matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错
//                                                            if (Utils.bzCheckPeptideWithAllChargeAtSameRT)
//                                                            {
//
//                                                            }
//                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);

                                                            if (ms2Int!=null) {
                                                                matchPep.dPrositRT = ms2Int.dPredictRT;
                                                            }
                                                            matchPep.calculateCombineScore();
                                                            matchPep.setDms2rt(arrMSTwoRts[iRTPos]);


                                                            matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);


                                                            if(Utils.bzCheckMS2TrailContainMS1Trail)
                                                            {
                                                                for (int iMassPos = 0; iMassPos < msonetrail.getArrMS1Mzs().length; iMassPos++) {
                                                                    for(int jr = 0;(window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos])!=null) && jr<window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).size();jr++)
                                                                    {
                                                                        if (msonetrail.getArrMS1Mzs()[iMassPos] == window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).get(jr).getMzApex()) {
                                                                            iMS1TrailInMS2TrailCount ++;
                                                                            break;
                                                                        }
                                                                    }

                                                                }
//                                                                matchPep.setStrMS1trailInMS2TrailCount(""+iMS1TrailInMS2TrailCount);


                                                            }
                                                            dPeptideMatchScoreSum += matchPep.dMatch;//将在时间内该Peptide所有的可以匹配的值全部加起来
                                                            if (matchPeptdideAll == null) {
                                                                matchPeptdideAll = matchPep;

                                                            } else {
                                                                matchPeptdideAll.combineResult(matchPep);
//                                                                matchPeptdideAll.combineResultSetRepeatCoverCount(matchPep);//update repeat and coverd
                                                            }
                                                            if (matchPep.dMatch > dPeptideMatchMax) {
                                                                dPeptideMatchMax = matchPep.dMatch;
                                                                matchPeptdideMax = matchPep;
                                                            }
                                                        }

////节省时间不能保存值 begin
//                                                    double score = 0.0;
//                                                    if (mapTriplePeptideRTScore[iPos]==null)
//                                                    {
////                                                        System.out.println("----"+arrMSTwoRts[iRTPos]);
//
//                                                        score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
//                                                        mapTriplePeptideRTScore[iPos] = new BitSet();
//                                                        mapTriplePeptideRTScore[iPos].set(iRTPos);
//                                                    }else
//                                                    {
//                                                        if(mapTriplePeptideRTScore[iPos].get(iRTPos)==false)
//                                                        {
//                                                            score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
////                                                            System.out.println("----~~~~~~~~~"+arrMSTwoRts[iRTPos]);
//
//                                                            mapTriplePeptideRTScore[iPos].set(iRTPos);
//
//                                                        }
//                                                    }
////节省时间不能保存值 end

//                                                    if(score>=0.4) {
//                                                        try {
//                                                            bufferedWriter.write(icount.incrementAndGet() + "\t");
//                                                            bufferedWriter.write(peps.get(iPos).id+ "\t" +peps.get(iPos).composition + "\t"+ peps.get(iPos).mass);
//                                                            bufferedWriter.write( "\t"+msonetrail.getMass());
//                                                            bufferedWriter.write("\t"+score);
//
//                                                            bufferedWriter.write("\t"+arrMSTwoRts[iRTPos]+"\n");
//                                                        } catch (IOException e) {
//                                                            e.printStackTrace();
//                                                        }
//
////                                                        System.out.println(icount.incrementAndGet() +"--------");
////                                                        System.out.println(peps.get(iPos).composition + peps.get(iPos).mass);
////                                                        System.out.println(msonetrail.getMass());
////                                                        System.out.println(score);
//
//                                                        System.out.println(arrMSTwoRts[iRTPos]);
//                                                    }
                                                    }
                                                    MatchPeptideWithMultiCharge mp = null;
//                                                    MatchPeptdide mp = null;
                                                    double dPm;
                                                    if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
                                                    {
                                                        mp = matchPeptdideMax;
                                                        dPm = dPeptideMatchMax;
                                                        if(mp!=null)
                                                        {
//                                                            mp.setRepeatCoverCount(matchPeptdideAll);
                                                            if(Utils.bzCheckMS2TrailContainMS1Trail) {
                                                                mp.setStrMS1trailInMS2TrailCount("" + iMS1TrailInMS2TrailCount);
                                                            }
                                                        }/*else
                                                        {
                                                            //TEST BUG
                                                            System.out.println(matchPeptdideAll);

                                                        }*/

                                                    } else {
                                                        mp = matchPeptdideAll;
                                                        dPm = dPeptideMatchScoreSum;
                                                    }


                                                    //peptide search mode: the more ion covered, the more score have
                                                    if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.CoverMax) {

                                                        if (mp != null) {
                                                            int length = mp.lBIonCovered.stream().distinct().toArray().length +
                                                                    mp.lYIonCovered.stream().distinct().toArray().length;
//                                                        mp.dMatch = mp.dMatch * length;
                                                            mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                                                        }

                                                    } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.Repeat) {
                                                        if (mp != null) {
                                                            int length = mp.iBRepeatCount + mp.iYRepeatCount;
                                                            mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                                                        }

                                                    } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.RepeatPlusCoverMax) {
                                                        if (mp != null) {
                                                            int length = mp.iBRepeatCount + mp.iYRepeatCount + mp.lBIonCovered.stream().distinct().toArray().length +
                                                                    mp.lYIonCovered.stream().distinct().toArray().length;

                                                            mp.dMatch = mp.dMatch * length;//    / mp.strPepcomposition.length();
                                                        }

                                                    }


                                                    //remain candidate list
                                                    int iCandidateListSize = Utils.iCandidateListSize;//Candidate list no
                                                    if (matchPeptdideAll != null) {
                                                        if (listMaxiumPeptide.size() == 0) {
                                                            listMaxiumPeptide.add(mp);
                                                        } else {
                                                            int iL = 0;
                                                            Boolean bzAdd = false;
                                                            for (; iL < iCandidateListSize && iL < listMaxiumPeptide.size(); iL++) {
                                                                if (listMaxiumPeptide.get(iL).dMatch < dPm) {
                                                                    listMaxiumPeptide.add(iL, mp);//strategy one

                                                                    bzAdd = true;
                                                                    break;

                                                                }
                                                            }
//                                                            if (!bzAdd && iL < iCandidateListSize - 1) {
                                                            if (!bzAdd && iL < iCandidateListSize ) {

                                                                listMaxiumPeptide.add(mp);

                                                            }
                                                            if (listMaxiumPeptide.size() > iCandidateListSize) {
                                                                listMaxiumPeptide.removeLast();
                                                            }


                                                        }
                                                    }

                                                }//end of mz of peptide
                                            }
                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);
                                        }
//                                        }


                                        //get all peps mass equals the msonetrail
//                                        peps.stream()
//                                                .filter(e->Math.abs(e.mass-msonetrail.getMass())<Utils.thresholdPPM*e.mass)
//                                                .forEach(
//                                                        pep->{
//                                                            //compare pep with all mstwotrail in rt window by using b y ions
//                                                            uniqueValues.stream().filter(v->Math.abs(y-v) <= Utils.thresholdRT).forEach(
//                                                                    twoRTtime->{
////                                                                         if(DbMatch.PepMatch(pep,mapRTMSTwoTrail.get(twoRTtime),Config.scoreMode)>0.0) {
////                                                                             System.out.println(icount.incrementAndGet() +"--------");
////                                                                             System.out.println(pep.composition + pep.mass);
////                                                                             System.out.println(msonetrail.getMass());
////                                                                             System.out.println(twoRTtime);
////                                                                         }
//
////                                                                        mapRTMSTwoTrail.get(twoRTtime).stream().forEach(
////                                                                                msTwoTrail -> {
////                                                                                    DbMatch.PepMatch(pep,msTwoTrail,Config.scoreMode);
////
////
////                                                                                }
////                                                                        );
//                                                                    }
//                                                            );
//
//
//                                                        }
//                                                );
//                                    }//end of the window with mz of ms1 feature
                                    }
                            );
                        }





//                        i.getAndIncrement();
//                        long time = System.currentTimeMillis() - start;
 //                       System.out.println("#########time:" + i + "/" + arrRts.length);

                        //                       System.out.println(time);
                    }
            );
            pb.stop();




//            //根据RT，和预测的RT的拟合值，找出在RT范围内的PEPTIDE，添加precursor到现有的msoneTRAIL中。
//            //设置这个precursor只检查这个PEPTIDE(以前是检查所有peptide，导致非常慢，可以将rt慢长一些[-6.5,+6.5]，加1-4个电荷)
//            //考虑不同的电荷检查raw data里面是否存在误差容限内的isotope是否存在，可以检查一到两个isotope：1电荷是+1，2电荷是+0.5，3电荷是+0.33，4电荷是+0.25
//            //get  RT distance between ms1 interval 得到一个rt间隔的所有PEPTIDE
//
////            int iAutoRTPosBegin = Arrays.binarySearch(dpepsAutoRT, y-RTWidnowDistance[0]);
////            iAutoRTPosBegin = Utils.getiPosBeginOfNearest(dpepsAutoRT.length, iAutoRTPosBegin);
////            if (iAutoRTPosBegin > 0)
////                iAutoRTPosBegin = iAutoRTPosBegin - 1;//move to less than position
//            for (int iRTPos = 0;
//                 iRTPos < dpepsAutoRT.length ;
//                 iRTPos++) {
//
//                String strPep = strpepsAutoRT.get(iRTPos);
//                double dmass = Peptide.CalculateMassByPeptideStr(strPep);
//                for (int icharge = 1;icharge<=Utils.iChargeRange;icharge++)
//                {
//                    double dmz = Utils.MassToMz(dmass,icharge);
//
//                    //直接匹配相对应的PEPTIDE
//                    MSTwoTrailSet window = iswc.FindWindowWithMZ(dmz, 1);
////                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
//                    // if mass matches, then pep match
//                    if (window != null) {
//                        double[] arrMSTwoRts = window.arrMSTwoRts;
//                        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;
//
//                        double rtbegin = dpepsAutoRT[iRTPos]-Utils.thresholdPredictRT;
//
//                        int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, rtbegin);
//                        iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
//                        if (iTwoRTPosBegin > 0)
//                            iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
//
//
//                        for (int iMS2RTPos = iTwoRTPosBegin;
//                             iMS2RTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < dpepsAutoRT[iRTPos] + Utils.thresholdPredictRT;
//                             iMS2RTPos++) {
//                            //只匹配一个PEPTIDE,但基于多个时间
//                            peps.get(iPos).GenerateIons();
//                            double[][] arrIntensity = null;
//                            MS2Intensity ms2Int = MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ());
//                            if (ms2Int!=null)
//                                arrIntensity = ms2Int.arrdIntensity;
//
//                            MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiCharge(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                    Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail
//                                    ,arrIntensity);
//                                                       /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
//                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
//                                                                msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
//                            if (matchPep != null) {
//                                matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错
////                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);
//
//                                if (ms2Int != null)
//                                    matchPep.dPrositRT = ms2Int.dPredictRT;
//                                matchPep.calculateCombineScore();
//                                matchPep.setDms2rt(arrMSTwoRts[iRTPos]);
//
//
//                                matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);
//                            }
//
//                        }
////                        MSOneTrail tmpPrecursor = new MSOneTrail(
////                                iPrecursorNO.getAndIncrement(), dmz, dpepsAutoRT[iRTPos],
////                                icharge, 0, 0,
////                                0, 0, 0,
////                                0, 0, 0,
////                                0, String.valueOf(dmz), String.valueOf(dpepsAutoRT[iRTPos]), "0"
////
////                        );
////                        tmpPrecursor.strGeneratePeptide = strPep;
////                        tmpPrecursor.bAddPrecursor = true;//set the precurosor bz
////                        lMsoneTrail.add(tmpPrecursor);
//                    }
//                }
//
//            }

//            for (Peptide pep: peps) {
//                // search for the window
//                IsolationWindow window = spec.FindWindow(pep.mass, 2);
//                // if mass matches, then pep match
//                if (window != null) {
//                    DbMatch res = DbMatch.PepMatch(pep, window, Config.scoreMode);
//                    if (res != null) {
//                        if (res.matchedRts.score > 0) {
//                            resLst.add(res);
//                        }
//                    }
//                }
//            }

//            List<MatchPeptdide> sorted = new LinkedList<>();

            List<MatchPeptideWithMultiCharge> sorted = new LinkedList<>();

//            List<MatchPeptdide> sorted = mapMSOneFeaturePeptide.entrySet().stream()
////                    .sorted(Comparator.comparing(e -> e.getValue().stream().map(MatchPeptdide::getdMatch).min(Comparator.reverseOrder()).orElse((double) 0)))
//                    //and also sort each group before collecting them in one list
//                    .flatMap(e -> e.getValue().stream().sorted(Comparator.comparing(MatchPeptdide::getdMatch))).collect(Collectors.toList());
//
//            for (MatchPeptdide mp:sorted)
//                {
//                    bufferedWriter.write(mp.toString()+"\n");
//
//                }
            //
            //
            Set<Long> keys = mapMSOneFeaturePeptide.keySet();

//            if (Utils.bzJustOneTrailCanMatchAPeptide) Utils.OutputMode = Utils.OutputModeEnum.OutputALL;//consider all ms1feature can match peptide, because those match score should be adjust

            for (Long k : keys) {

                List<MatchPeptideWithMultiCharge> lmp = mapMSOneFeaturePeptide.get(k);
//                List<MatchPeptdide> lmp = mapMSOneFeaturePeptide.get(k);

                if (lmp.size() > 0)
                {
                    if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                    {
                        sorted.addAll(lmp);

                    }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
                    {
                        sorted.add(lmp.get(0));
                    }
                    else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
                    {
                        if (lmp.size() > 1) {
                            sorted.add(lmp.get(1));
                        }
                    }
                }

//                sorted.addAll(lmp);

            }



           /* if(Utils.bzCalculateThePearsonCorresionWithPredictMSMS)
            {
                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致

                *//*FileReader freader;
                BufferedReader br;
                String dataFile = Config.spectrumFolder +Utils.strpDeepOutfile;//pdeep3output.txt";
                String line;
                try {
                    freader = new FileReader(dataFile);
                    br = new BufferedReader(freader);
                    while ((line = br.readLine()) != null && !(line.startsWith(">peptide|")));//定位到当前位置

                } catch (FileNotFoundException noFile) {
                    throw new FileNotFoundException();
                }*//*


                for(int iS = 0 ;iS<sorted.size();iS++) {
//                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS_JustAll(
                            MapMS2Intensity.get(
                                    sorted.get(iS).strPepcomposition+sorted.get(iS).msOneTrail.getZ()).arrdIntensity);
                    if(iS%10000==0)
                        System.out.println(iS + "---------iS---------" + sorted.get(iS).getdMatch() + sorted.get(iS).strDecoy);

                   *//* //open predict MSMS results file
                    double[][] arrPredictMSMS = new double[sorted.get(iS).pep.composition.length() - 1][4];//length -1 with b b+2 y y+2

                    String[] strPredictPepInfo = line.split("\\|");

                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序

                    while (!strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        if ((line = br.readLine()) != null && (line.startsWith(">peptide|"))) {
                            strPredictPepInfo = line.split("\\|");
                        }
                    }
                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序，过滤掉不匹配的数据

                    if (strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        String strValue = "";
                        while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) {
                            strValue = strValue + line;
                        }
                        String[] arrstrValue = strValue.replaceAll("\\[ ", "")
                                .replaceAll("\\[", "")
                                .replaceAll("\\]", "")
                                .replaceAll("  ", " ")
                                .split(" ");
                        for (int iP = 0; iP < strPredictPepInfo[1].length() - 1; iP++) {
                            for (int jP = 0; jP < 4; jP++) {
                                if(arrstrValue[iP * 4 + jP].isEmpty())
                                {
                                    System.out.println(sorted.get(iS).pep.composition);//test
                                }else
                                {
                                    arrPredictMSMS[iP][jP] = Double.parseDouble(arrstrValue[iP * 4 + jP]);

                                }
                            }
                        }
                    }
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(arrPredictMSMS);*//*
                }
               // br.close();
               // freader.close();
                //predict msms results end
            }
            System.out.println("---------CalculateThePearsonCorresionWithPredictMSMS Finish---------");
*/
            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long,Integer> mapIDCount = new HashMap();

            if (Utils.bzJustOneTrailCanMatchAPeptide) //need to adjust the match score
            {

                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致


                for(int iS = 0 ;iS<20000 && iS<sorted.size();iS++) {
                        int jS = iS + 1;
                    if (sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) {
                                l = mid + 1;
                            } else {
                                h = mid;
                            }
                        }

                        if (l - iS > 1 && sorted.get(iS).getdMatch() < sorted.get(iS + 1).getdMatch()) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS + "----jS--------------" + (jS < sorted.size() ? sorted.get(jS).getdMatch() : "last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        } else {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }
                    } else {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }
                }

                 /*   if(sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) l = mid + 1;
                            else h = mid;
                        }


//                        for(;jS<sorted.size();jS++)
//                        {
//                            if(sorted.get(iS).getdMatch() >= sorted.get(jS).getdMatch()) break;
//                        }
                        if(jS-iS>1) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS+"----jS--------------"+(jS<sorted.size()?sorted.get(jS).getdMatch():"last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        }else
                        {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }

                    }else
                    {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }

                }*/
/*                for (MatchPeptdide mp : sorted) {
                    mp.adjustScore();
                }*/


            }
            Collections.sort(sorted, Collections.reverseOrder());

            //output multiple file, concise file, only peptidefile
            FileWriter concisefileWriter;
            BufferedWriter concisebufferedWriter = null;
            FileWriter peptidesForpDeep3Writer = null;
            BufferedWriter peptidesForpDeep3bufferedWriter = null;
            FileWriter peptidesForpPrositWriter = null;
            BufferedWriter peptidesForpPrositbufferedWriter = null;
            if(Utils.bzOutputConciseAndPeptideFiles)
            {
                 concisefileWriter = new FileWriter(Config.spectrumFolder + "concisefile_0502.csv");
                 concisebufferedWriter = new BufferedWriter(concisefileWriter);
                peptidesForpDeep3Writer = new FileWriter(Config.spectrumFolder + "peptidesForpDeep3_0502.csv");
                peptidesForpDeep3bufferedWriter = new BufferedWriter(peptidesForpDeep3Writer);
                peptidesForpPrositWriter = new FileWriter(Config.spectrumFolder + "peptidesForpProsit_0502.csv");
                peptidesForpPrositbufferedWriter = new BufferedWriter(peptidesForpPrositWriter);
                concisebufferedWriter.write("msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tgetQuality_score\tmsonemz\tmsonecharge\tBIon\tYIon"+"\n");
                peptidesForpDeep3bufferedWriter.write("peptide\tmodinfo\tcharge"+"\n");
                peptidesForpPrositbufferedWriter.write("modified_sequence\tcollision_energy\tprecursor_charge\tfragmentation\n");

            }
            System.out.println("---------OutputALL Begin---------");

            for (MatchPeptideWithMultiCharge mp : sorted) {
                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
//                bufferedWriter.write(k+ "\t");
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");

//                    if(Utils.bzOutputConciseAndPeptideFiles)
//                    {
//
//                        concisebufferedWriter.write(mp.getConciseString()+"\n");
//                        peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
//                        peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");
//
//                    }
                }else
                {
//                    System.out.println(mp.strPepcomposition);
                    bufferedWriter.write(mp.toString() + "\n");


                }
                if(Utils.bzOutputConciseAndPeptideFiles)
                {

                    concisebufferedWriter.write(mp.getConciseString()+"\n");
                    peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
                    peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");

                }

                iCountForFDR++;

            }
            System.out.println("---------OutputALL Finish---------");

            if(Utils.bzOutputConciseAndPeptideFiles) {
                concisebufferedWriter.flush();
                peptidesForpDeep3bufferedWriter.flush();
                peptidesForpPrositbufferedWriter.flush();
                concisebufferedWriter.close();
                peptidesForpDeep3bufferedWriter.close();
                peptidesForpPrositbufferedWriter.close();
            }
        }
//        System.out.println("The number of matched peps is: " + );


        bufferedWriter.flush();
        bufferedWriter.close();

        // sort result
//        Collections.sort(resLst, Collections.reverseOrder());
//        System.out.println("The number of matched peps is: " + resLst.size());
//        bufferedWriter.flush();
//        bufferedWriter.close();

//        for (DbMatch match : resLst) {
//            // TODO Remove all System.out.println calls for production
//            System.out.println(match.id);
//            System.out.println(match.composition);
//            System.out.println(match.matchedRts.score);
//
//            // Output the information to a file
//            bufferedWriter.write(match.id + "\n");
//            bufferedWriter.write(match.composition + "\n");
//            bufferedWriter.write(Double.toString(match.matchedRts.score) + "\n");
//
//            if (Config.mode == Enums.RunMode.DEBUG) {
//                System.out.println(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")) + "\n");
//            } else {
//                System.out.println(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")) + "\n");
//            }
//        }
    }


    public void match_TrailWithMultiRT(String mgfInfile, String fastaInfile, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();

        // read fasta
        ArrayList<Genome> genomes = FastaFile.ReadFile(fastaInfile);
//        ArrayList<Genome> genomes = FastaFile.ReadFileFormPeptide(fastaInfile);


        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top
        List<Peptide> peps = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
        double[] pepsMass = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());

//        FileWriter fileWriterPeptideInfo = new FileWriter("peptideInfo.csv");
//        BufferedWriter bufferedWriterPeptideInfo = new BufferedWriter(fileWriterPeptideInfo);
//        bufferedWriterPeptideInfo.write( "peptide\tmass\tdMutationRate"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
//                + "\n");
//        for (Peptide pep:peps)
//        {
//            bufferedWriterPeptideInfo.write(pep.composition+'\t'+pep.mass+'\t'+pep.dMutationRate+'\n');
//        }
//
//        bufferedWriterPeptideInfo.flush();
//        bufferedWriterPeptideInfo.close();


        ArrayList<DbMatch> resLst = new ArrayList<>();




//        Set<Peptide>
        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tmsoneMZ\tmsoneZ\tmsonePeakAreaLocalRank" +
                "\tms2rt\tpeptideLength\tMS1PeptideMassError\tdPeptideMutationRate\tibyComplemetaryCount\tdCosWithPredictMSMS\tdMS1WithMS2TrailCoscinSimilarity\tMS2TrailCoscinSimilarity\tBIon\tYIon");
        for(int iw =0;iw<Utils.iMS2Charge;iw++) {
            bufferedWriter.write("\tdlogBionPeakAreaSum_C"+(iw+1)+"\tdlogYionPeakAreaSum_C"+(iw+1)+"\tdBionMassErrorSum_C"+(iw+1)+"\tdYionMassErrorSum_C"+(iw+1)+"\tiBionPeakHalfSum_C"+(iw+1)+"\tiYionPeakHalfSum_C"+(iw+1)+"\tbionMatchCount_C"+(iw+1)+"\tyionMatchCount" +
                    "_C"+(iw+1)+"\tdBionRetentionTimeErrorSum_C"+(iw+1)+"\tdYionRetentionTimeErrorSum_C"+(iw+1)+"\tiBConsective_C"+(iw+1)+"\tiYConsective_C"+(iw+1)+"\tdBionCosWithPredictMSMS_C"+(iw+1)+"\tdYionCosWithPredictMSMS_C"+(iw+1));
        }
        bufferedWriter.write( "\tarrMatchedBion\tarrMatchedYion\tadjustBYIon\tSoredIndex\twindowSize\tdBionCos60WithPredictMSMS_C1\tdCosWithPredictMSMSMatchedPredict\tprositRT\tcountSize"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
                + "\n");



//        Map<Integer,Map<Integer,Double>> mapTriplePeptideRTScore = new TreeMap<>();

//        BitSet mapTriplePeptideRTScore = new BitSet(peps.size());
//        BitSet[] mapTriplePeptideRTScore = new BitSet[peps.size()];
//        Map<Pair<Integer,Integer>,Double> mapTriplePeptideRTScore = new HashMap<>();
//        for (MSTwoTrail msTwoTrail : MSTWOts.arrMSTwoTrail) {
//            double[] rts = msTwoTrail.getRts();
//
//            for(double drt:rts) {
//                if (uniqueValues.add(drt))
//                {
//                    List<MSTwoTrail> listMSTwoTrail = new ArrayList<>();
//                    listMSTwoTrail.add(msTwoTrail);
//
//                    mapRTMSTwoTrail.put(drt,listMSTwoTrail);
//                }
//                else
//                {
//                    mapRTMSTwoTrail.get(drt).add(msTwoTrail);
//                }
//            }
//        }

        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            spec.readFromeTrailFile(file);
            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();

//            Map<Double, List<MSOneTrail>> mapRTMSoneTrail =  spec.arrMSOneTrail.stream().
//                    collect(Collectors.groupingBy(x-> x.getRt()));

            TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));





            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();


            //GET MS1 feature trials from RT, which are adopted to remove tails from ms2 tails, because there are precusors that not be fragmented.

//            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges1.out",mapRTMSoneTrail);
            IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges.out",mapRTMSoneTrail);


//            Map<String, MS2Intensity> MapMS2Intensity = Main.getMapMS2Intensity();
            Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit();


            //set ms1feature peakarea rank
            int ims1 = 1;
            for(MSOneTrail ms1Trail:spec.arrMSOneTrail)
            {
                long iPeakAreaRank = 0;
                MSTwoTrailSet window = iswc.FindWindowWithMZ(ms1Trail.getMz(), 1);
                //get ms1feature RT distance between Math.abs(y-v) <= Utils.thresholdRT
                int iMS1RTPosBegin = Arrays.binarySearch(arrRts, ms1Trail.getRt() - Utils.thresholdMS1RT);
                iMS1RTPosBegin = Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);
                if (iMS1RTPosBegin > 0) {
                    iMS1RTPosBegin = iMS1RTPosBegin - 1;//move to less than position
                }

                for (int iRTPos = iMS1RTPosBegin;
                     iRTPos < arrRts.length && arrRts[iRTPos] < ms1Trail.getRt() + Utils.thresholdMS1RT;
                     iRTPos++) {
                    if(window!=null) {
                        iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
                                x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh && x.getQuantification_peaks_area() >= ms1Trail.getQuantification_peaks_area()
                        ).count();
                    }
                }
                if(iPeakAreaRank>0)
                {
                    ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20 / iPeakAreaRank), 0));
                }else
                {
                    ms1Trail.setdPeakAreaLocalRankInTimeSpan( 0 );
                }


/*                if(window!=null) {
                    iPeakAreaRank = spec.arrMSOneTrail.stream().filter((x) ->
                            x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh && x.getQuantification_peaks_area() >= ms1Trail.getQuantification_peaks_area()
                                    && x.getRt() >= ms1Trail.getRt() - Utils.thresholdMS1RT && x.getRt() <= ms1Trail.getRt() + Utils.thresholdMS1RT
                    ).count();
                }*/

//                System.out.println("ims1:"+ ims1++  +"----"+spec.arrMSOneTrail.size());

            }

//        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(
//                "isolationWindowRanges.txt"));
//        IsolationWindowCollection iswc = (IsolationWindowCollection) ois.readObject();

//        List<double[]> larrMSTwoRts = new ArrayList<>(iswTest.windows.size());
//        List<Map<Double, List<MSTwoTrail>>> lmapRTMStwoTrail = new ArrayList<>(iswTest.windows.size());

            // TODO: READ IN SPECTRUM RANKING FILE (Added by Caroline) Implement Rank Here


            // for each isolation window collection, match once
            //Step 1: get a msonetrail
            //Step 2: sort with RT
            //Step 3:for each msonetrail get the window
            //      check all peptide with same mass  (mz +- 30 ppm)
            //         for each mstwoTrail in window (mz in window && RT +- 5s )
            //             compare mstwotrail b y ions with peptide b y ions (mz +- 30 ppm)
            //         add peptide to msonetrail as candidate peptide with confident score
//        MSTwoTrailSet MSTWOts = new MSTwoTrailSet();
//        MSTWOts.readFromeTrailFile(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window0_ms2_trails.tsv");

//        Map<Double,List<MSTwoTrail> >mapRTMSTwoTrail = new HashMap<Double, List<MSTwoTrail>>();
            TreeSet<Double> uniqueValues = new TreeSet<>();


            //////////////////////////////////
            //for a window,initial rts and rts corresponding trails
//        double[] arrMSTwoRts = MSTWOts.arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
//                .distinct()
//                .sorted()
//                .toArray();
////        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().
////                collect(Collectors.groupingBy(x-> x.getRtApex()));
//
//
//        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
            //////////////////////////////////////

            Map<Long, LinkedList<MatchPeptideWithMultiCharge>> mapMSOneFeaturePeptide = new HashMap<>();

//            Map<Long, LinkedList<MatchPeptdide>> mapMSOneFeaturePeptide = new HashMap<>();

            AtomicInteger i = new AtomicInteger();
            Arrays.stream(arrRts).forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();

                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);
                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {


//                                        if (msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass >= pepsMass[0]) //filter m1 trail which less than the least peptide mass

                                        //begin of window with mz of m1 feature
////                                        MSTwoTrailSet window = iswc.FindWindow(msonetrail.getMass(), 2);
//                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                        // if mass matches, then pep match
//                                        if (window != null) {
//                                            double[] arrMSTwoRts = window.arrMSTwoRts;
//                                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;

                                        if (msonetrail.getMass() - Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {



                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window

                                            LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide = new LinkedList<MatchPeptideWithMultiCharge>();

//                                            LinkedList<MatchPeptdide> listMaxiumPeptide = new LinkedList<MatchPeptdide>();


                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况
//                                            while (iPosBegin > 0) {
//                                                iPosBegin = iPosBegin - 1;//move to less than position
//                                                if(iPosBegin > 0 && pepsMass[iPosBegin]>pepsMass[iPosBegin-1]) break;//PEPTIDE的质量可能有相等的情况，这是peptide有序情况
//                                            }

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {



                                                //begin of window with mz of ms1 feature
                                                MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
                                                // if mass matches, then pep match
                                                if (window != null) {
                                                    double[] arrMSTwoRts = window.arrMSTwoRts;
                                                    Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;


                                                    double dPeptideMatchScoreSum = 0.0;//用于记录每个peptide在整个保留时间的得分
                                                    double dPeptideMatchMax = -1000000.0;

                                                    int iMS1TrailInMS2TrailCount = 0;

                                                    MatchPeptideWithMultiCharge matchPeptdideMax = null;
                                                    MatchPeptideWithMultiCharge matchPeptdideAll = null;

                                                    // recode each b y ion for best match ion
                                                    MSTwoTrail[] BionMatchTrail = new MSTwoTrail[peps.get(iPos).composition.length()*Utils.iMS2Charge];
                                                    MSTwoTrail[] YionMatchTrail = new MSTwoTrail[peps.get(iPos).composition.length()*Utils.iMS2Charge];



                                                    //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
                                                    int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getRt() - Utils.thresholdRT);
                                                    iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
                                                    if (iTwoRTPosBegin > 0) {
                                                        iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
                                                    }



                                                    for (int iRTPos = iTwoRTPosBegin;
                                                         iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getRt() + Utils.thresholdRT;
                                                         iRTPos++) {

                                                        //TEST_DEBUG
                                                       /* if(DebugRTTime(msonetrail.getRt(),peps.get(iPos).composition))
                                                        {
                                                            System.out.println(msonetrail.getRt()+"  "+peps.get(iPos).composition);
                                                            try {
                                                                FileWriter fileWriterTest = new FileWriter(psmOutfile+"html/"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass()+".html");
                                                                BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterTest);
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Title+"_"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass());
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi1Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));

                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi2Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart2Ms2trailData);
                                                                int lenth = mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size();
                                                                double[][] ddData= new double[lenth][3];
//                                                                bufferedWriterTest.write(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass()+"\n");
                                                                int iL =0;
                                                                for(MSTwoTrail mstwot:mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]))
                                                                {
                                                                    ddData[iL][0]=mstwot.getMzApex();
                                                                    ddData[iL][1]=mstwot.getPeakArea();
                                                                    ddData[iL][2]=iL+1;
                                                                    iL++;


                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().mapToDouble(MSTwoTrail::getMzApex).toArray())+"\n");

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddData)+"\n");




                                                                double  ponDA = 1.00727647;

                                                                double[][] ddBData= new double[peps.get(iPos).b_ions.length][3];


                                                                for(int ibs=0;ibs<peps.get(iPos).b_ions.length;ibs++)
                                                                {
                                                                    ddBData[ibs][0]=peps.get(iPos).b_ions[ibs];
                                                                    ddBData[ibs][1]=500;
                                                                    ddBData[ibs][2]=ibs+1;


                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3BData);


                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");
                                                                double [] b_Charge = peps.get(iPos).b_ions.clone();
                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)
                                                                {
                                                                    ddBData[ibs][0]=(b_Charge[ibs]+ponDA)/2;//已经加了一个
                                                                    ddBData[ibs][1]=400;
                                                                    ddBData[ibs][2]=ibs+1;
//                                                                    b_Charge[ibs] = (b_Charge[ibs]+2*ponDA)/2;
                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3B2ChargeData);

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");

//                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)//3charge
//                                                                {
//                                                                    ddBData[ibs][0]=(ddBData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddBData[ibs][1]=300;
//                                                                    ddBData[ibs][2]=ibs+1;
////                                                                    b_Charge[ibs] = (b_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");


                                                                double[][] ddYData= new double[peps.get(iPos).y_ions.length][3];


                                                                for(int ibs=0;ibs<peps.get(iPos).y_ions.length;ibs++)
                                                                {
                                                                    ddYData[ibs][0]=peps.get(iPos).y_ions[ibs];
                                                                    ddYData[ibs][1]=500;
//                                                                    ddYData[ibs][2]=ibs+1;
                                                                    ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3YData);

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                                                                double [] y_Charge = peps.get(iPos).y_ions.clone();
                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
                                                                {
                                                                    ddYData[ibs][0]=(y_Charge[ibs]+ponDA)/2;//已经加了一个
                                                                    ddYData[ibs][1]=400;
//                                                                    ddYData[ibs][2]=ibs+1;
                                                                    ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;


//                                                                    y_Charge[ibs]=(y_Charge[ibs]+2*ponDA)/2;
                                                                }
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3Y2ChargeData);

                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(y_Charge)+"\n");
//                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
//                                                                {
//                                                                    ddYData[ibs][0]=(ddYData[ibs][0]*2+ponDA)/3;//已经加了两个
//                                                                    ddYData[ibs][1]=300;
//                                                                    ddYData[ibs][2]=ibs+1;
////                                                                    y_Charge[ibs]=(y_Charge[ibs]*2+ponDA)/3;
//                                                                }
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart4End);

                                                                bufferedWriterTest.flush();
                                                                bufferedWriterTest.close();

                                                            } catch (IOException e) {
                                                                e.printStackTrace();
                                                            }
                                                        }*/

                                                        ////为获得最大值，或者是加和值，每一次重复计算
                                                        double score = 0.0;

//                                                       MatchPeptdide matchPep = DbMatch.PepMatch(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]), Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(), msonetrail.getId());

                                                        //    System.out.println(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass());

                                                        peps.get(iPos).GenerateIons();
                                                        double[][] arrIntensity = null;
                                                        MS2Intensity ms2Int = MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ());
                                                        if (ms2Int!=null) {
                                                            arrIntensity = ms2Int.arrdIntensity;
                                                        }

                                                        MatchPeptideWithMultiCharge matchPep = DbMatch.PepMatchIonCountCrossRTMultiCharge(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail
                                                                ,arrIntensity);
                                                       /* MatchPeptdide matchPep = DbMatch.PepMatchIonCountCrossRT(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                                                                Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(),
                                                                msonetrail.getId(),msonetrail.getQuality_score(),BionMatchTrail,YionMatchTrail);*/
                                                        if (matchPep != null) {
                                                            matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错
//                                                            matchPep.calculatePearsonCorrelationWithPredictMSMS_JustAll(MapMS2Intensity.get(peps.get(iPos).composition+msonetrail.getZ()).arrdIntensity);

                                                            if (ms2Int!=null) {
                                                                matchPep.dPrositRT = ms2Int.dPredictRT;
                                                            }
                                                            matchPep.calculateCombineScore();
                                                            matchPep.setDms2rt(arrMSTwoRts[iRTPos]);


                                                            matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);


                                                            if(Utils.bzCheckMS2TrailContainMS1Trail)
                                                            {
                                                                for (int iMassPos = 0; iMassPos < msonetrail.getArrMS1Mzs().length; iMassPos++) {
                                                                    for(int jr = 0;(window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos])!=null) && jr<window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).size();jr++)
                                                                    {
                                                                        if (msonetrail.getArrMS1Mzs()[iMassPos] == window.mapRTMStwoTrailFromRemoveMS1Trail.get(arrMSTwoRts[iRTPos]).get(jr).getMzApex()) {
                                                                            iMS1TrailInMS2TrailCount ++;
                                                                            break;
                                                                        }
                                                                    }

                                                                }
//                                                                matchPep.setStrMS1trailInMS2TrailCount(""+iMS1TrailInMS2TrailCount);


                                                            }
                                                            dPeptideMatchScoreSum += matchPep.dMatch;//将在时间内该Peptide所有的可以匹配的值全部加起来
                                                            if (matchPeptdideAll == null) {
                                                                matchPeptdideAll = matchPep;

                                                            } else {
                                                                matchPeptdideAll.combineResult(matchPep);
//                                                                matchPeptdideAll.combineResultSetRepeatCoverCount(matchPep);//update repeat and coverd
                                                            }
                                                            if (matchPep.dMatch > dPeptideMatchMax) {
                                                                dPeptideMatchMax = matchPep.dMatch;
                                                                matchPeptdideMax = matchPep;
                                                            }
                                                        }

////节省时间不能保存值 begin
//                                                    double score = 0.0;
//                                                    if (mapTriplePeptideRTScore[iPos]==null)
//                                                    {
////                                                        System.out.println("----"+arrMSTwoRts[iRTPos]);
//
//                                                        score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
//                                                        mapTriplePeptideRTScore[iPos] = new BitSet();
//                                                        mapTriplePeptideRTScore[iPos].set(iRTPos);
//                                                    }else
//                                                    {
//                                                        if(mapTriplePeptideRTScore[iPos].get(iRTPos)==false)
//                                                        {
//                                                            score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
////                                                            System.out.println("----~~~~~~~~~"+arrMSTwoRts[iRTPos]);
//
//                                                            mapTriplePeptideRTScore[iPos].set(iRTPos);
//
//                                                        }
//                                                    }
////节省时间不能保存值 end

//                                                    if(score>=0.4) {
//                                                        try {
//                                                            bufferedWriter.write(icount.incrementAndGet() + "\t");
//                                                            bufferedWriter.write(peps.get(iPos).id+ "\t" +peps.get(iPos).composition + "\t"+ peps.get(iPos).mass);
//                                                            bufferedWriter.write( "\t"+msonetrail.getMass());
//                                                            bufferedWriter.write("\t"+score);
//
//                                                            bufferedWriter.write("\t"+arrMSTwoRts[iRTPos]+"\n");
//                                                        } catch (IOException e) {
//                                                            e.printStackTrace();
//                                                        }
//
////                                                        System.out.println(icount.incrementAndGet() +"--------");
////                                                        System.out.println(peps.get(iPos).composition + peps.get(iPos).mass);
////                                                        System.out.println(msonetrail.getMass());
////                                                        System.out.println(score);
//
//                                                        System.out.println(arrMSTwoRts[iRTPos]);
//                                                    }
                                                    }
                                                    MatchPeptideWithMultiCharge mp = null;
//                                                    MatchPeptdide mp = null;
                                                    double dPm;
                                                    if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
                                                    {
                                                        mp = matchPeptdideMax;
                                                        dPm = dPeptideMatchMax;
                                                        if(mp!=null)
                                                        {
//                                                            mp.setRepeatCoverCount(matchPeptdideAll);
                                                            if(Utils.bzCheckMS2TrailContainMS1Trail) {
                                                                mp.setStrMS1trailInMS2TrailCount("" + iMS1TrailInMS2TrailCount);
                                                            }
                                                        }/*else
                                                        {
                                                            //TEST BUG
                                                            System.out.println(matchPeptdideAll);

                                                        }*/

                                                    } else {
                                                        mp = matchPeptdideAll;
                                                        dPm = dPeptideMatchScoreSum;
                                                    }


                                                    //peptide search mode: the more ion covered, the more score have
                                                    if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.CoverMax) {

                                                        if (mp != null) {
                                                            int length = mp.lBIonCovered.stream().distinct().toArray().length +
                                                                    mp.lYIonCovered.stream().distinct().toArray().length;
//                                                        mp.dMatch = mp.dMatch * length;
                                                            mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                                                        }

                                                    } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.Repeat) {
                                                        if (mp != null) {
                                                            int length = mp.iBRepeatCount + mp.iYRepeatCount;
                                                            mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                                                        }

                                                    } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.RepeatPlusCoverMax) {
                                                        if (mp != null) {
                                                            int length = mp.iBRepeatCount + mp.iYRepeatCount + mp.lBIonCovered.stream().distinct().toArray().length +
                                                                    mp.lYIonCovered.stream().distinct().toArray().length;
                                                            ;
                                                            mp.dMatch = mp.dMatch * length;//    / mp.strPepcomposition.length();
                                                        }

                                                    }


                                                    //remain candidate list
                                                    int iCandidateListSize = Utils.iCandidateListSize;//Candidate list no
                                                    if (matchPeptdideAll != null) {
                                                        if (listMaxiumPeptide.size() == 0) {
                                                            listMaxiumPeptide.add(mp);
                                                        } else {
                                                            int iL = 0;
                                                            Boolean bzAdd = false;
                                                            for (; iL < iCandidateListSize && iL < listMaxiumPeptide.size(); iL++) {
                                                                if (listMaxiumPeptide.get(iL).dMatch < dPm) {
                                                                    listMaxiumPeptide.add(iL, mp);//strategy one

                                                                    bzAdd = true;
                                                                    break;

                                                                }
                                                            }
//                                                            if (!bzAdd && iL < iCandidateListSize - 1) {
                                                            if (!bzAdd && iL < iCandidateListSize ) {

                                                                listMaxiumPeptide.add(mp);

                                                            }
                                                            if (listMaxiumPeptide.size() > iCandidateListSize) {
                                                                listMaxiumPeptide.removeLast();
                                                            }


                                                        }
                                                    }

                                                }//end of mz of peptide
                                            }
                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);
                                        }
//                                        }


                                        //get all peps mass equals the msonetrail
//                                        peps.stream()
//                                                .filter(e->Math.abs(e.mass-msonetrail.getMass())<Utils.thresholdPPM*e.mass)
//                                                .forEach(
//                                                        pep->{
//                                                            //compare pep with all mstwotrail in rt window by using b y ions
//                                                            uniqueValues.stream().filter(v->Math.abs(y-v) <= Utils.thresholdRT).forEach(
//                                                                    twoRTtime->{
////                                                                         if(DbMatch.PepMatch(pep,mapRTMSTwoTrail.get(twoRTtime),Config.scoreMode)>0.0) {
////                                                                             System.out.println(icount.incrementAndGet() +"--------");
////                                                                             System.out.println(pep.composition + pep.mass);
////                                                                             System.out.println(msonetrail.getMass());
////                                                                             System.out.println(twoRTtime);
////                                                                         }
//
////                                                                        mapRTMSTwoTrail.get(twoRTtime).stream().forEach(
////                                                                                msTwoTrail -> {
////                                                                                    DbMatch.PepMatch(pep,msTwoTrail,Config.scoreMode);
////
////
////                                                                                }
////                                                                        );
//                                                                    }
//                                                            );
//
//
//                                                        }
//                                                );
//                                    }//end of the window with mz of ms1 feature
                                    }
                            );
                        }
                        i.getAndIncrement();
                        long time = System.currentTimeMillis() - start;
                        System.out.println("#########time:" + i + "/" + arrRts.length);

                        //                       System.out.println(time);
                    }
            );


//            for (Peptide pep: peps) {
//                // search for the window
//                IsolationWindow window = spec.FindWindow(pep.mass, 2);
//                // if mass matches, then pep match
//                if (window != null) {
//                    DbMatch res = DbMatch.PepMatch(pep, window, Config.scoreMode);
//                    if (res != null) {
//                        if (res.matchedRts.score > 0) {
//                            resLst.add(res);
//                        }
//                    }
//                }
//            }

//            List<MatchPeptdide> sorted = new LinkedList<>();

            List<MatchPeptideWithMultiCharge> sorted = new LinkedList<>();

//            List<MatchPeptdide> sorted = mapMSOneFeaturePeptide.entrySet().stream()
////                    .sorted(Comparator.comparing(e -> e.getValue().stream().map(MatchPeptdide::getdMatch).min(Comparator.reverseOrder()).orElse((double) 0)))
//                    //and also sort each group before collecting them in one list
//                    .flatMap(e -> e.getValue().stream().sorted(Comparator.comparing(MatchPeptdide::getdMatch))).collect(Collectors.toList());
//
//            for (MatchPeptdide mp:sorted)
//                {
//                    bufferedWriter.write(mp.toString()+"\n");
//
//                }
            //
            //
            Set<Long> keys = mapMSOneFeaturePeptide.keySet();

//            if (Utils.bzJustOneTrailCanMatchAPeptide) Utils.OutputMode = Utils.OutputModeEnum.OutputALL;//consider all ms1feature can match peptide, because those match score should be adjust

            for (Long k : keys) {

                List<MatchPeptideWithMultiCharge> lmp = mapMSOneFeaturePeptide.get(k);
//                List<MatchPeptdide> lmp = mapMSOneFeaturePeptide.get(k);

                if (lmp.size() > 0)
                {
                    if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                    {
                        sorted.addAll(lmp);

                    }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
                    {
                        sorted.add(lmp.get(0));
                    }
                    else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
                    {
                        if (lmp.size() > 1) {
                            sorted.add(lmp.get(1));
                        }
                    }
                }

//                sorted.addAll(lmp);

            }



           /* if(Utils.bzCalculateThePearsonCorresionWithPredictMSMS)
            {
                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致

                *//*FileReader freader;
                BufferedReader br;
                String dataFile = Config.spectrumFolder +Utils.strpDeepOutfile;//pdeep3output.txt";
                String line;
                try {
                    freader = new FileReader(dataFile);
                    br = new BufferedReader(freader);
                    while ((line = br.readLine()) != null && !(line.startsWith(">peptide|")));//定位到当前位置

                } catch (FileNotFoundException noFile) {
                    throw new FileNotFoundException();
                }*//*


                for(int iS = 0 ;iS<sorted.size();iS++) {
//                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS_JustAll(
                            MapMS2Intensity.get(
                                    sorted.get(iS).strPepcomposition+sorted.get(iS).msOneTrail.getZ()).arrdIntensity);
                    if(iS%10000==0)
                        System.out.println(iS + "---------iS---------" + sorted.get(iS).getdMatch() + sorted.get(iS).strDecoy);

                   *//* //open predict MSMS results file
                    double[][] arrPredictMSMS = new double[sorted.get(iS).pep.composition.length() - 1][4];//length -1 with b b+2 y y+2

                    String[] strPredictPepInfo = line.split("\\|");

                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序

                    while (!strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        if ((line = br.readLine()) != null && (line.startsWith(">peptide|"))) {
                            strPredictPepInfo = line.split("\\|");
                        }
                    }
                    //若预测的序列是ms2 spectrum的子集，且是一样的顺序，过滤掉不匹配的数据

                    if (strPredictPepInfo[1].equals(sorted.get(iS).pep.composition)) {
                        String strValue = "";
                        while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) {
                            strValue = strValue + line;
                        }
                        String[] arrstrValue = strValue.replaceAll("\\[ ", "")
                                .replaceAll("\\[", "")
                                .replaceAll("\\]", "")
                                .replaceAll("  ", " ")
                                .split(" ");
                        for (int iP = 0; iP < strPredictPepInfo[1].length() - 1; iP++) {
                            for (int jP = 0; jP < 4; jP++) {
                                if(arrstrValue[iP * 4 + jP].isEmpty())
                                {
                                    System.out.println(sorted.get(iS).pep.composition);//test
                                }else
                                {
                                    arrPredictMSMS[iP][jP] = Double.parseDouble(arrstrValue[iP * 4 + jP]);

                                }
                            }
                        }
                    }
                    sorted.get(iS).calculatePearsonCorrelationWithPredictMSMS(arrPredictMSMS);*//*
                }
               // br.close();
               // freader.close();
                //predict msms results end
            }
            System.out.println("---------CalculateThePearsonCorresionWithPredictMSMS Finish---------");
*/
            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long,Integer> mapIDCount = new HashMap();

            if (Utils.bzJustOneTrailCanMatchAPeptide) //need to adjust the match score
            {

                Collections.sort(sorted, Collections.reverseOrder());//保证输出序列与已预测的MSMS数据序列一致


                for(int iS = 0 ;iS<20000 && iS<sorted.size();iS++) {
                    int jS = iS + 1;
                    if (sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) {
                                l = mid + 1;
                            } else {
                                h = mid;
                            }
                        }

                        if (l - iS > 1 && sorted.get(iS).getdMatch() < sorted.get(iS + 1).getdMatch()) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS + "----jS--------------" + (jS < sorted.size() ? sorted.get(jS).getdMatch() : "last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        } else {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }
                    } else {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }
                }

                 /*   if(sorted.get(iS).adjustScore())//调整了得分需要向后移动
                    {
                        int l = jS;
                        int h = sorted.size();
                        while (l < h) {
                            int mid = (l + h) / 2;
                            if (sorted.get(mid).getdMatch() <= sorted.get(iS).getdMatch()) l = mid + 1;
                            else h = mid;
                        }


//                        for(;jS<sorted.size();jS++)
//                        {
//                            if(sorted.get(iS).getdMatch() >= sorted.get(jS).getdMatch()) break;
//                        }
                        if(jS-iS>1) {

                            sorted.add(l, sorted.get(iS));
                            sorted.remove(sorted.get(iS));


                            System.out.println(jS+"----jS--------------"+(jS<sorted.size()?sorted.get(jS).getdMatch():"last"));

//                            Collections.rotate(sorted.subList(iS, jS), -1);
                            iS--;
                        }else
                        {
                            sorted.get(iS).setMS2TrailSelected();//不需要调整位置就设置ms2TRAIL被占用

                        }

                    }else
                    {
                        sorted.get(iS).setMS2TrailSelected();//不需要调整就设置ms2TRAIL被占用
                    }

                }*/
/*                for (MatchPeptdide mp : sorted) {
                    mp.adjustScore();
                }*/


            }
            Collections.sort(sorted, Collections.reverseOrder());

            //output multiple file, concise file, only peptidefile
            FileWriter concisefileWriter;
            BufferedWriter concisebufferedWriter = null;
            FileWriter peptidesForpDeep3Writer = null;
            BufferedWriter peptidesForpDeep3bufferedWriter = null;
            FileWriter peptidesForpPrositWriter = null;
            BufferedWriter peptidesForpPrositbufferedWriter = null;
            if(Utils.bzOutputConciseAndPeptideFiles)
            {
                concisefileWriter = new FileWriter(Config.spectrumFolder + "concisefile_0502.csv");
                concisebufferedWriter = new BufferedWriter(concisefileWriter);
                peptidesForpDeep3Writer = new FileWriter(Config.spectrumFolder + "peptidesForpDeep3_0502.csv");
                peptidesForpDeep3bufferedWriter = new BufferedWriter(peptidesForpDeep3Writer);
                peptidesForpPrositWriter = new FileWriter(Config.spectrumFolder + "peptidesForpProsit_0502.csv");
                peptidesForpPrositbufferedWriter = new BufferedWriter(peptidesForpPrositWriter);
                concisebufferedWriter.write("msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tgetQuality_score\tmsonemz\tmsonecharge\tBIon\tYIon"+"\n");
                peptidesForpDeep3bufferedWriter.write("peptide\tmodinfo\tcharge"+"\n");
                peptidesForpPrositbufferedWriter.write("modified_sequence\tcollision_energy\tprecursor_charge\tfragmentation\n");

            }
            System.out.println("---------OutputALL Begin---------");

            for (MatchPeptideWithMultiCharge mp : sorted) {
                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
//                bufferedWriter.write(k+ "\t");
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");

//                    if(Utils.bzOutputConciseAndPeptideFiles)
//                    {
//
//                        concisebufferedWriter.write(mp.getConciseString()+"\n");
//                        peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
//                        peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");
//
//                    }
                }else
                {
//                    System.out.println(mp.strPepcomposition);
                    bufferedWriter.write(mp.toString() + "\n");


                }
                if(Utils.bzOutputConciseAndPeptideFiles)
                {

                    concisebufferedWriter.write(mp.getConciseString()+"\n");
                    peptidesForpDeep3bufferedWriter.write(mp.pep.composition+"\t\t"+mp.msOneTrail.getZ()+"\n");
                    peptidesForpPrositbufferedWriter.write(mp.pep.composition+"\t27.5\t"+mp.msOneTrail.getZ()+"\tHCD\n");

                }

                iCountForFDR++;

            }
            System.out.println("---------OutputALL Finish---------");

            if(Utils.bzOutputConciseAndPeptideFiles) {
                concisebufferedWriter.flush();
                peptidesForpDeep3bufferedWriter.flush();
                peptidesForpPrositbufferedWriter.flush();
                concisebufferedWriter.close();
                peptidesForpDeep3bufferedWriter.close();
                peptidesForpPrositbufferedWriter.close();
            }
        }
//        System.out.println("The number of matched peps is: " + );


        bufferedWriter.flush();
        bufferedWriter.close();

        // sort result
//        Collections.sort(resLst, Collections.reverseOrder());
//        System.out.println("The number of matched peps is: " + resLst.size());
//        bufferedWriter.flush();
//        bufferedWriter.close();

//        for (DbMatch match : resLst) {
//            // TODO Remove all System.out.println calls for production
//            System.out.println(match.id);
//            System.out.println(match.composition);
//            System.out.println(match.matchedRts.score);
//
//            // Output the information to a file
//            bufferedWriter.write(match.id + "\n");
//            bufferedWriter.write(match.composition + "\n");
//            bufferedWriter.write(Double.toString(match.matchedRts.score) + "\n");
//
//            if (Config.mode == Enums.RunMode.DEBUG) {
//                System.out.println(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")) + "\n");
//            } else {
//                System.out.println(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")) + "\n");
//            }
//        }
    }

/*
    public void match(String mgfInfile, String fastaInfile, String psmOutfile) throws Throwable {
        // TODO: Possibly refactor mgfInfile to be a string array of different filenames

        String[] files = {mgfInfile};
        AtomicInteger icount = new AtomicInteger();

        // read fasta
        ArrayList<Genome> genomes = FastaFile.ReadFile(fastaInfile);
//        ArrayList<Genome> genomes = FastaFile.ReadFileFormPeptide(fastaInfile);


        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top
        List<Peptide> peps = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
        double[] pepsMass = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= Utils.thresholdPeptideSizeMin).distinct().sorted().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());

        ArrayList<DbMatch> resLst = new ArrayList<>();


//        //Load the windows
        IsolationWindowCollection iswc = new IsolationWindowCollection(Config.spectrumFolder + "isolationWindowRanges.out");

//        ObjectInputStream ois = new ObjectInputStream(new FileInputStream(
//                "isolationWindowRanges.txt"));
//        IsolationWindowCollection iswc = (IsolationWindowCollection) ois.readObject();

//        List<double[]> larrMSTwoRts = new ArrayList<>(iswTest.windows.size());
//        List<Map<Double, List<MSTwoTrail>>> lmapRTMStwoTrail = new ArrayList<>(iswTest.windows.size());

        // TODO: READ IN SPECTRUM RANKING FILE (Added by Caroline) Implement Rank Here


        // for each isolation window collection, match once
        //Step 1: get a msonetrail
        //Step 2: sort with RT
        //Step 3:for each msonetrail get the window
        //      check all peptide with same mass  (mz +- 30 ppm)
        //         for each mstwoTrail in window (mz in window && RT +- 5s )
        //             compare mstwotrail b y ions with peptide b y ions (mz +- 30 ppm)
        //         add peptide to msonetrail as candidate peptide with confident score
//        MSTwoTrailSet MSTWOts = new MSTwoTrailSet();
//        MSTWOts.readFromeTrailFile(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window0_ms2_trails.tsv");

//        Map<Double,List<MSTwoTrail> >mapRTMSTwoTrail = new HashMap<Double, List<MSTwoTrail>>();
        TreeSet<Double> uniqueValues = new TreeSet<>();


        //////////////////////////////////
        //for a window,initial rts and rts corresponding trails
//        double[] arrMSTwoRts = MSTWOts.arrMSTwoTrail.stream().mapToDouble(MSTwoTrail::getRtApex)
//                .distinct()
//                .sorted()
//                .toArray();
////        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().
////                collect(Collectors.groupingBy(x-> x.getRtApex()));
//
//
//        Map<Double, List<MSTwoTrail>> mapRTMStwoTrail =  MSTWOts.arrMSTwoTrail.stream().sorted()
//                .collect(Collectors.groupingBy(MSTwoTrail::getRtApex, TreeMap::new, Collectors.toList()));
        //////////////////////////////////////

//        Set<Peptide>
        FileWriter fileWriter = new FileWriter(psmOutfile);
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tms2rt\tBIon\tYIon\t" +
                "dlogBionPeakAreaSum\tdlogYionPeakAreaSum\tdBionMassErrorSum\tdYionMassErrorSum\tiBionPeakHalfSum\tiYionPeakHalfSum\tbionMatchCount\tyionMatchCount\t" +
                "peptideLength\tMS1PeptideMassError\tiBConsective\tiYConsective\tiBRepeatCount\tiYRepeatCount\tiBCoverCount\tiYCoverCount\tbStrCovered\tyStrCovered"
                + "\n");



//        Map<Integer,Map<Integer,Double>> mapTriplePeptideRTScore = new TreeMap<>();

//        BitSet mapTriplePeptideRTScore = new BitSet(peps.size());
//        BitSet[] mapTriplePeptideRTScore = new BitSet[peps.size()];
//        Map<Pair<Integer,Integer>,Double> mapTriplePeptideRTScore = new HashMap<>();
//        for (MSTwoTrail msTwoTrail : MSTWOts.arrMSTwoTrail) {
//            double[] rts = msTwoTrail.getRts();
//
//            for(double drt:rts) {
//                if (uniqueValues.add(drt))
//                {
//                    List<MSTwoTrail> listMSTwoTrail = new ArrayList<>();
//                    listMSTwoTrail.add(msTwoTrail);
//
//                    mapRTMSTwoTrail.put(drt,listMSTwoTrail);
//                }
//                else
//                {
//                    mapRTMSTwoTrail.get(drt).add(msTwoTrail);
//                }
//            }
//        }

        for (String file : files) {

            MSOneTrailSet spec = new MSOneTrailSet();
            spec.readFromeTrailFile(file);
            double[] arrRts = spec.arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                    .distinct()
                    .sorted()
                    .toArray();

//            Map<Double, List<MSOneTrail>> mapRTMSoneTrail =  spec.arrMSOneTrail.stream().
//                    collect(Collectors.groupingBy(x-> x.getRt()));

            Map<Double, List<MSOneTrail>> mapRTMSoneTrail = spec.arrMSOneTrail.stream().sorted()
                    .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));


            double RTWidnowDistance[] = IntStream.range(0, arrRts.length - 1)
                    .mapToDouble(i -> {
                        return (arrRts[i + 1] - arrRts[i]);
                    })
                    .toArray();


            Map<Long, LinkedList<MatchPeptdide>> mapMSOneFeaturePeptide = new HashMap<>();

            AtomicInteger i = new AtomicInteger();
            Arrays.stream(arrRts).forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();

                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);
                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {


//                                        if (msonetrail.getMass() - Utils.thresholdPPM * pepsMass[0] - Utils.H2OMass >= pepsMass[0]) //filter m1 trail which less than the least peptide mass

                                        //begin of window with mz of m1 feature
////                                        MSTwoTrailSet window = iswc.FindWindow(msonetrail.getMass(), 2);
//                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                        // if mass matches, then pep match
//                                        if (window != null) {
//                                            double[] arrMSTwoRts = window.arrMSTwoRts;
//                                            Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;

                                        if (msonetrail.getMass() - Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {



                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window


                                            LinkedList<MatchPeptdide> listMaxiumPeptide = new LinkedList<MatchPeptdide>();


                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况
//                                            while (iPosBegin > 0) {
//                                                iPosBegin = iPosBegin - 1;//move to less than position
//                                                if(iPosBegin > 0 && pepsMass[iPosBegin]>pepsMass[iPosBegin-1]) break;//PEPTIDE的质量可能有相等的情况，这是peptide有序情况
//                                            }

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {



                                                //begin of window with mz of ms1 feature
                                                MSTwoTrailSet window = iswc.FindWindowWithMZ(msonetrail.getMz(), 1);
//                                                MSTwoTrailSet window = iswc.FindWindow(pepsMass[iPos], 2);
                                                // if mass matches, then pep match
                                                if (window != null) {
                                                    double[] arrMSTwoRts = window.arrMSTwoRts;
                                                    Map<Double, List<MSTwoTrail>> mapRTMStwoTrail = window.mapRTMStwoTrail;


                                                    double dPeptideMatchScoreSum = 0.0;//用于记录每个peptide在整个保留时间的得分
                                                    double dPeptideMatchMax = -1.0;
                                                    MatchPeptdide matchPeptdideMax = null;
                                                    MatchPeptdide matchPeptdideAll = null;



                                                    //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
                                                    int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getRt() - Utils.thresholdRT);
                                                    iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
                                                    if (iTwoRTPosBegin > 0)
                                                        iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position



                                                    for (int iRTPos = iTwoRTPosBegin;
                                                         iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getRt() + Utils.thresholdRT;
                                                         iRTPos++) {

//                                                        //TEST_DEBUG
//                                                        if(DebugRTTime(msonetrail.getRt(),peps.get(iPos).composition))
//                                                        {
//                                                            System.out.println(msonetrail.getRt()+"  "+peps.get(iPos).composition);
//                                                            try {
//                                                                FileWriter fileWriterTest = new FileWriter(psmOutfile+"html/"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass()+".html");
//                                                                BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterTest);
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Title+"_"+peps.get(iPos).composition+"_MS1RT:"+msonetrail.getRt()+"_MS2RT:"+arrMSTwoRts[iRTPos]+"_"+msonetrail.getMass());
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi1Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));
//
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart1Axi2Max+(Math.max(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).get(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size()-1).getMzApex(),peps.get(iPos).mass)+50));
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart2Ms2trailData);
//                                                                int lenth = mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).size();
//                                                                double[][] ddData= new double[lenth][3];
////                                                                bufferedWriterTest.write(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass()+"\n");
//                                                                int iL =0;
//                                                                for(MSTwoTrail mstwot:mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]))
//                                                                {
//                                                                    ddData[iL][0]=mstwot.getMzApex();
//                                                                    ddData[iL][1]=mstwot.getPeakArea();
//                                                                    ddData[iL][2]=iL+1;
//                                                                    iL++;
//
//
//                                                                }
////                                                                bufferedWriterTest.write(JSONWriter.valueToString(mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]).stream().mapToDouble(MSTwoTrail::getMzApex).toArray())+"\n");
//
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddData)+"\n");
//
//
//
//
//                                                                double  ponDA = 1.00727647;
//
//                                                                double[][] ddBData= new double[peps.get(iPos).b_ions.length][3];
//
//
//                                                                for(int ibs=0;ibs<peps.get(iPos).b_ions.length;ibs++)
//                                                                {
//                                                                    ddBData[ibs][0]=peps.get(iPos).b_ions[ibs];
//                                                                    ddBData[ibs][1]=500;
//                                                                    ddBData[ibs][2]=ibs+1;
//
//
//                                                                }
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3BData);
//
//
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");
//                                                                double [] b_Charge = peps.get(iPos).b_ions.clone();
//                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)
//                                                                {
//                                                                    ddBData[ibs][0]=(b_Charge[ibs]+ponDA)/2;//已经加了一个
//                                                                    ddBData[ibs][1]=400;
//                                                                    ddBData[ibs][2]=ibs+1;
////                                                                    b_Charge[ibs] = (b_Charge[ibs]+2*ponDA)/2;
//                                                                }
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3B2ChargeData);
//
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");
//
////                                                                for(int ibs=0;ibs<b_Charge.length;ibs++)//3charge
////                                                                {
////                                                                    ddBData[ibs][0]=(ddBData[ibs][0]*2+ponDA)/3;//已经加了两个
////                                                                    ddBData[ibs][1]=300;
////                                                                    ddBData[ibs][2]=ibs+1;
//////                                                                    b_Charge[ibs] = (b_Charge[ibs]*2+ponDA)/3;
////                                                                }
////                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddBData)+"\n");
//
//
//                                                                double[][] ddYData= new double[peps.get(iPos).y_ions.length][3];
//
//
//                                                                for(int ibs=0;ibs<peps.get(iPos).y_ions.length;ibs++)
//                                                                {
//                                                                    ddYData[ibs][0]=peps.get(iPos).y_ions[ibs];
//                                                                    ddYData[ibs][1]=500;
////                                                                    ddYData[ibs][2]=ibs+1;
//                                                                    ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;
//
//
//                                                                }
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3YData);
//
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
//                                                                double [] y_Charge = peps.get(iPos).y_ions.clone();
//                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
//                                                                {
//                                                                    ddYData[ibs][0]=(y_Charge[ibs]+ponDA)/2;//已经加了一个
//                                                                    ddYData[ibs][1]=400;
////                                                                    ddYData[ibs][2]=ibs+1;
//                                                                    ddYData[ibs][2]=peps.get(iPos).y_ions.length-ibs;
//
//
////                                                                    y_Charge[ibs]=(y_Charge[ibs]+2*ponDA)/2;
//                                                                }
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart3Y2ChargeData);
//
//                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
////                                                                bufferedWriterTest.write(JSONWriter.valueToString(y_Charge)+"\n");
////                                                                for(int ibs=0;ibs<y_Charge.length;ibs++)
////                                                                {
////                                                                    ddYData[ibs][0]=(ddYData[ibs][0]*2+ponDA)/3;//已经加了两个
////                                                                    ddYData[ibs][1]=300;
////                                                                    ddYData[ibs][2]=ibs+1;
//////                                                                    y_Charge[ibs]=(y_Charge[ibs]*2+ponDA)/3;
////                                                                }
////                                                                bufferedWriterTest.write(JSONWriter.valueToString(ddYData)+"\n");
//                                                                bufferedWriterTest.write(Utils.strHtmlTemplatePart4End);
//
//                                                                bufferedWriterTest.flush();
//                                                                bufferedWriterTest.close();
//
//                                                            } catch (IOException e) {
//                                                                e.printStackTrace();
//                                                            }
//                                                        }

                                                        ////为获得最大值，或者是加和值，每一次重复计算
                                                        double score = 0.0;

//                                                       MatchPeptdide matchPep = DbMatch.PepMatch(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]), Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(), msonetrail.getId());

                                                    //    System.out.println(peps.get(iPos).composition+"\t"+arrMSTwoRts[iRTPos]+"\t"+msonetrail.getRt()+"\t"+msonetrail.getMass());

                                                        peps.get(iPos).GenerateIons();
                                                        MatchPeptdide matchPep = DbMatch.PepMatchIonCount(peps.get(iPos), mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]), Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail.getRt(), msonetrail.getMass(), msonetrail.getId(),msonetrail.getQuality_score());
                                                        if (matchPep != null) {
                                                            matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错
                                                            matchPep.calculateCombineScore();
                                                            matchPep.setDms2rt(arrMSTwoRts[iRTPos]);
                                                            dPeptideMatchScoreSum += matchPep.dMatch;//将在时间内该Peptide所有的可以匹配的值全部加起来
                                                            if (matchPeptdideAll == null) {
                                                                matchPeptdideAll = matchPep;

                                                            } else {
                                                                matchPeptdideAll.combineResult(matchPep);
                                                                matchPeptdideAll.combineResultSetRepeatCoverCount(matchPep);//update repeat and coverd
                                                            }
                                                            if (matchPep.dMatch > dPeptideMatchMax) {
                                                                dPeptideMatchMax = matchPep.dMatch;
                                                                matchPeptdideMax = matchPep;
                                                            }
                                                        }

////节省时间不能保存值 begin
//                                                    double score = 0.0;
//                                                    if (mapTriplePeptideRTScore[iPos]==null)
//                                                    {
////                                                        System.out.println("----"+arrMSTwoRts[iRTPos]);
//
//                                                        score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
//                                                        mapTriplePeptideRTScore[iPos] = new BitSet();
//                                                        mapTriplePeptideRTScore[iPos].set(iRTPos);
//                                                    }else
//                                                    {
//                                                        if(mapTriplePeptideRTScore[iPos].get(iRTPos)==false)
//                                                        {
//                                                            score = DbMatch.PepMatch(peps.get(iPos),mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]),Config.scoreMode,bufferedWriter,icount.getAndIncrement(),msonetrail.getRt(),msonetrail.getMass());
////                                                            System.out.println("----~~~~~~~~~"+arrMSTwoRts[iRTPos]);
//
//                                                            mapTriplePeptideRTScore[iPos].set(iRTPos);
//
//                                                        }
//                                                    }
////节省时间不能保存值 end

//                                                    if(score>=0.4) {
//                                                        try {
//                                                            bufferedWriter.write(icount.incrementAndGet() + "\t");
//                                                            bufferedWriter.write(peps.get(iPos).id+ "\t" +peps.get(iPos).composition + "\t"+ peps.get(iPos).mass);
//                                                            bufferedWriter.write( "\t"+msonetrail.getMass());
//                                                            bufferedWriter.write("\t"+score);
//
//                                                            bufferedWriter.write("\t"+arrMSTwoRts[iRTPos]+"\n");
//                                                        } catch (IOException e) {
//                                                            e.printStackTrace();
//                                                        }
//
////                                                        System.out.println(icount.incrementAndGet() +"--------");
////                                                        System.out.println(peps.get(iPos).composition + peps.get(iPos).mass);
////                                                        System.out.println(msonetrail.getMass());
////                                                        System.out.println(score);
//
//                                                        System.out.println(arrMSTwoRts[iRTPos]);
//                                                    }
                                                    }

                                                    MatchPeptdide mp = null;
                                                    double dPm;
                                                    if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
                                                    {
                                                        mp = matchPeptdideMax;
                                                        dPm = dPeptideMatchMax;
                                                        if(mp!=null)
                                                            mp.setRepeatCoverCount(matchPeptdideAll);

                                                    } else {
                                                        mp = matchPeptdideAll;
                                                        dPm = dPeptideMatchScoreSum;
                                                    }


                                                    //peptide search mode: the more ion covered, the more score have
                                                    if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.CoverMax) {

                                                        if (mp != null) {
                                                            int length = mp.lBIonCovered.stream().distinct().toArray().length +
                                                                    mp.lYIonCovered.stream().distinct().toArray().length;
//                                                        mp.dMatch = mp.dMatch * length;
                                                            mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                                                        }

                                                    } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.Repeat) {
                                                        if (mp != null) {
                                                            int length = mp.iBRepeatCount + mp.iYRepeatCount;
                                                            mp.dMatch = mp.dMatch * length;//     / mp.strPepcomposition.length();
                                                        }

                                                    } else if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.RepeatPlusCoverMax) {
                                                        if (mp != null) {
                                                            int length = mp.iBRepeatCount + mp.iYRepeatCount + mp.lBIonCovered.stream().distinct().toArray().length +
                                                                    mp.lYIonCovered.stream().distinct().toArray().length;
                                                            ;
                                                            mp.dMatch = mp.dMatch * length;//    / mp.strPepcomposition.length();
                                                        }

                                                    }


                                                    //remain candidate list
                                                    int iCandidateListSize = Utils.iCandidateListSize;//Candidate list no
                                                    if (matchPeptdideAll != null) {
                                                        if (listMaxiumPeptide.size() == 0) {
                                                            listMaxiumPeptide.add(mp);
                                                        } else {
                                                            int iL = 0;
                                                            Boolean bzAdd = false;
                                                            for (; iL < iCandidateListSize && iL < listMaxiumPeptide.size(); iL++) {
                                                                if (listMaxiumPeptide.get(iL).dMatch < dPm) {
                                                                    listMaxiumPeptide.add(iL, mp);//strategy one

                                                                    bzAdd = true;
                                                                    break;

                                                                }
                                                            }
                                                            if (!bzAdd && iL < iCandidateListSize - 1) {
                                                                listMaxiumPeptide.add(mp);

                                                            }
                                                            if (listMaxiumPeptide.size() > iCandidateListSize)
                                                                listMaxiumPeptide.removeLast();


                                                        }
                                                    }

                                                }//end of mz of peptide
                                            }
                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);
                                        }
//                                        }


                                        //get all peps mass equals the msonetrail
//                                        peps.stream()
//                                                .filter(e->Math.abs(e.mass-msonetrail.getMass())<Utils.thresholdPPM*e.mass)
//                                                .forEach(
//                                                        pep->{
//                                                            //compare pep with all mstwotrail in rt window by using b y ions
//                                                            uniqueValues.stream().filter(v->Math.abs(y-v) <= Utils.thresholdRT).forEach(
//                                                                    twoRTtime->{
////                                                                         if(DbMatch.PepMatch(pep,mapRTMSTwoTrail.get(twoRTtime),Config.scoreMode)>0.0) {
////                                                                             System.out.println(icount.incrementAndGet() +"--------");
////                                                                             System.out.println(pep.composition + pep.mass);
////                                                                             System.out.println(msonetrail.getMass());
////                                                                             System.out.println(twoRTtime);
////                                                                         }
//
////                                                                        mapRTMSTwoTrail.get(twoRTtime).stream().forEach(
////                                                                                msTwoTrail -> {
////                                                                                    DbMatch.PepMatch(pep,msTwoTrail,Config.scoreMode);
////
////
////                                                                                }
////                                                                        );
//                                                                    }
//                                                            );
//
//
//                                                        }
//                                                );
//                                    }//end of the window with mz of ms1 feature
                                    }
                            );
                        }
                        i.getAndIncrement();
                        long time = System.currentTimeMillis() - start;
                        System.out.println("#########time:" + i + "/" + arrRts.length);

 //                       System.out.println(time);
                    }
            );


//            for (Peptide pep: peps) {
//                // search for the window
//                IsolationWindow window = spec.FindWindow(pep.mass, 2);
//                // if mass matches, then pep match
//                if (window != null) {
//                    DbMatch res = DbMatch.PepMatch(pep, window, Config.scoreMode);
//                    if (res != null) {
//                        if (res.matchedRts.score > 0) {
//                            resLst.add(res);
//                        }
//                    }
//                }
//            }


            List<MatchPeptdide> sorted = new LinkedList<>();

//            List<MatchPeptdide> sorted = mapMSOneFeaturePeptide.entrySet().stream()
////                    .sorted(Comparator.comparing(e -> e.getValue().stream().map(MatchPeptdide::getdMatch).min(Comparator.reverseOrder()).orElse((double) 0)))
//                    //and also sort each group before collecting them in one list
//                    .flatMap(e -> e.getValue().stream().sorted(Comparator.comparing(MatchPeptdide::getdMatch))).collect(Collectors.toList());
//
//            for (MatchPeptdide mp:sorted)
//                {
//                    bufferedWriter.write(mp.toString()+"\n");
//
//                }
            //
            //
            Set<Long> keys = mapMSOneFeaturePeptide.keySet();
            for (Long k : keys) {
                List<MatchPeptdide> lmp = mapMSOneFeaturePeptide.get(k);

                if (lmp.size() > 0)
                {
                    if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                    {
                        sorted.addAll(lmp);

                    }else if(Utils.OutputMode == Utils.OutputModeEnum.OutputMax)
                    {
                        sorted.add(lmp.get(0));
                    }
                    else if(Utils.OutputMode == Utils.OutputModeEnum.OutputSecondMax)
                    {
                        if (lmp.size() > 1) {
                            sorted.add(lmp.get(1));
                        }
                    }
                }

//                sorted.addAll(lmp);

            }

            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long,Integer> mapIDCount = new HashMap();
            Collections.sort(sorted, Collections.reverseOrder());
            for (MatchPeptdide mp : sorted) {
                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
//                bufferedWriter.write(k+ "\t");
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if(Utils.OutputMode == Utils.OutputModeEnum.OutputALL)
                {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");
                }else
                {
                    bufferedWriter.write(mp.toString() + "\n");
                }




                iCountForFDR++;

            }
        }
//        System.out.println("The number of matched peps is: " + );


        bufferedWriter.flush();
        bufferedWriter.close();

        // sort result
//        Collections.sort(resLst, Collections.reverseOrder());
//        System.out.println("The number of matched peps is: " + resLst.size());
//        bufferedWriter.flush();
//        bufferedWriter.close();

//        for (DbMatch match : resLst) {
//            // TODO Remove all System.out.println calls for production
//            System.out.println(match.id);
//            System.out.println(match.composition);
//            System.out.println(match.matchedRts.score);
//
//            // Output the information to a file
//            bufferedWriter.write(match.id + "\n");
//            bufferedWriter.write(match.composition + "\n");
//            bufferedWriter.write(Double.toString(match.matchedRts.score) + "\n");
//
//            if (Config.mode == Enums.RunMode.DEBUG) {
//                System.out.println(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(Object::toString)
//                        .collect(Collectors.joining("\n")) + "\n");
//            } else {
//                System.out.println(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")));
//                bufferedWriter.write(match.matchedRts.peaks.stream().map(x -> x.rt.toString())
//                        .collect(Collectors.joining(" ")) + "\n");
//            }
//        }
    }

*/


    private boolean DebugRTTime(double rt,String strPeptide) {
        /*
        81.1575
        86.271
        93.97083333
        50.5985
        105.1156667
        112.9468333
         */

        double[] testRT = {
                88.9853333333333,
                86.271,
                93.97083333,
                50.5985,
                105.1156667,
                112.9468333,
                88.99816667,
                87.291,
                123.875,
                37.356,
                123.3128333,
                99.57133333,
                73.86283333,
                115.3201667,
                115.1311667,
                120.6605,
                46.053,
                87.30666667
        };
        String[] testStrPeptide = {
                "DFVAEPMGEKPVGSLAGIGEVLGK",
                "ISLGLPVGAVINCADNTGAK",
                "ALAPTWEQLALGLEHSETVK",
                "TYADYESVNECMEGVCK",
                "EPLFGISTGNLITGLAAGAK",
                "VGLTSEILNSFEHEFLSK",
                "ILTVEDHYYEGGIGEAVSSAVVGEPGITVTHLAVNR",
                "HIADLAGNSEVILPVPAFNVINGGSHAGNK",
                "TLSTIATSTDAASVVHSTDLVVEAIVENLK",
                "NMGGPYGGGNYGPGGSGGSGGYGGR",
                "FQSSAVMALQEASEAYLVGLFEDTNLCAIHAK",
                "LEDLSESIVNDFAYMK",
                "LTPLILKPFGNSISR",
                "APVVLLSEPACAHALEALATLLRPR",
                "TNLEFLQEQFNSIAAHVLHCTDSGFGAR",
                "SNVKPNSGELDPLYVVEVLLR",
                "LVNHFVEEFK",
                "ESVLTATSILNNPIVK",
        };
        int i=0;
        for(double  tRT:testRT)
        {
            if(rt>tRT-0.5 && rt<tRT+0.5 && strPeptide.equals(testStrPeptide[i])) {
                return true;
            }
            i++;
        }
        return false;

    }


}
