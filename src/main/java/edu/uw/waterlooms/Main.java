package edu.uw.waterlooms;

import edu.uw.waterlooms.match.*;
import edu.uw.waterlooms.match.DbMatch;
import edu.uw.waterlooms.match.FastaFile;
import edu.uw.waterlooms.match.Genome;
import edu.uw.waterlooms.match.Peptide;
import edu.uw.waterlooms.match.RunFilter;
import edu.uw.waterlooms.ms1.*;
import edu.uw.waterlooms.entity.*;

import java.io.*;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import edu.uw.waterlooms.msutil.MZXMLReader;
import edu.uw.waterlooms.msutil.OpenMzxml;
import edu.uw.waterlooms.peptideMatch.*;
import edu.uw.waterlooms.peptideMatch.Config;
import edu.uw.waterlooms.peptideMatch.Utils;
import edu.uw.waterlooms.service.ParameterService;
import me.tongfei.progressbar.ProgressBar;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.util.Pair;
import org.json.JSONArray;
import org.json.JSONObject;
import org.json.JSONWriter;

/**
 * WaterlooMS Main Class
 *
 * @author Shaokai Wang, Jia Wu, Xiangyuan Zeng, Ziwei Han
 */
public class Main {


    public static void main(String[] args) throws IOException {

        // Parse environment variable set in Dockerfile to validate if the executable is running within
        // the container
        boolean containerExecution = Boolean.parseBoolean(System.getenv("RUN_WITHIN_CONTAINER"));
        ArrayList<String> processBuilderCommand = new ArrayList<String>();
        if (!containerExecution) {
            processBuilderCommand.add("docker-compose");
            processBuilderCommand.add("exec");
            processBuilderCommand.add("-T");
            processBuilderCommand.add("waterlooms");
        }
        processBuilderCommand.add("python");

        Path path = FileSystems.getDefault().getPath("").toAbsolutePath();
        String DOCKER_WORKING_DIR = "/data/";
        String LOCAL_WORKING_DIR = containerExecution ? "/data/" : path + "/data/";

        ArgumentParser parser =
                ArgumentParsers.newFor("java -jar waterlooms.jar")
                        .build()
                        .defaultHelp(true)
                        .description("DIA Analysis.");
        parser
                .addArgument("-mzXML")
                .metavar("FILE")
                .type(Arguments.fileType().acceptSystemIn().verifyExists())
                .help("mzXML File");
        parser.addArgument("-parameters").type(String.class).help("JSON String Dict of Parameters");
        parser
                .addArgument("-fasta")
                .metavar("FILE")
                .type(Arguments.fileType().acceptSystemIn().verifyExists())
                .help("Fasta Peptide Database File. Must contain decoys.");
        parser
                .addArgument("-outputDir")
                .type(String.class)
                .help("Directory for outputting the result CSV.");


        parser
                .addArgument("-autort")
                .type(String.class)
                .help("Autort result file.");


        parser
                .addArgument("-prositDir")
                .type(String.class)
                .help("Prosit directory should include prosit result files.");
        // Parse Arguments
        Namespace ns = null;
        try {
            ns = parser.parseArgs(args);
        } catch (ArgumentParserException $exception) {
            parser.handleError($exception);
            System.exit(1);
        }

        String mzXMLInFile = ns.getString("mzXML");
        String fastaInFile = ns.getString("fasta");
        String outputDir = ns.getString("outputDir");
        String prositDir = ns.getString("prositDir");
        String autortFile = ns.getString("autort");
        String serializedParameters = ns.getString("parameters");

        // Configure the file for local debugging
        String rawFileName = "toy.mzXML";


        if ((mzXMLInFile != null && !mzXMLInFile.isEmpty())) {
            // TODO: DIA-WEBAPP submits mzXMLInFile as a path/file.mzXML
            // TODO: Need to strip this and set rawFileName
            rawFileName = mzXMLInFile;
        }


        // Set defaults from main/resources/ folder if arguments not specified (development purposes)
        String mzXMLFile =
                (mzXMLInFile != null && !mzXMLInFile.isEmpty())
                        ? mzXMLInFile
                        : LOCAL_WORKING_DIR + rawFileName;


        //follow mzxmlFile path set the docker_working_path
        String filepath = mzXMLFile.substring(0, mzXMLFile.lastIndexOf("/"));
        DOCKER_WORKING_DIR = filepath + "/";
        LOCAL_WORKING_DIR = filepath + "/";
        System.out.println(DOCKER_WORKING_DIR);


//
//        // test another docker and run r code
//        ArrayList<String> processBuilderRCommandTest = new ArrayList<String>();
////        processBuilderRCommand.add("docker");
////        processBuilderRCommand.add("exec");
////        processBuilderRCommand.add("-t");
////        processBuilderRCommand.add("rmodel_dp");//need change
//
//        processBuilderRCommandTest.add("Rscript");
//
//        System.out.println("Running model.r ...");
//        ProcessBuilder processBuilderRLTest = new ProcessBuilder();
//        ArrayList<String> rlCommandtest = new ArrayList<>(processBuilderRCommandTest);
//        rlCommandtest.add("/waterlooms/src/main/R/model.R" );
//        rlCommandtest.add("HumanR01_1125");
//        rlCommandtest.add("Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_psm.csv");
//        rlCommandtest.add(autortFile);//("autoRTHumanR01_1t3d.csv");
//        rlCommandtest.add(autortFile+"_minfit_maxfit_diff_1t3d.txt");//("humanR01_RTminfit_maxfit_diff_1t3d.txt");
//        rlCommandtest.add(DOCKER_WORKING_DIR);
//        processBuilderRLTest.command(rlCommandtest);
//        try {
//            // TODO: Suppress output if necessary for python
//            Process process = processBuilderRLTest.inheritIO().start();
//            int exitCode = process.waitFor();
//            System.out.println("\nmodel.r exited with code : " + exitCode);
//        } catch (InterruptedException $e) {
//            $e.printStackTrace();
//        }
//
////        docker exec -t lucid_hawking /usr/local/lib/R/bin/Rscript /home/rstudio/data/model.R
//
//        System.exit(1);


        String fastaFile =
                (fastaInFile != null && !fastaInFile.isEmpty())
                        ? fastaInFile
//            : LOCAL_WORKING_DIR + "database_withdecoy.fasta";
                        : LOCAL_WORKING_DIR + "homo_sapiens_withdecoy.fasta";
        String outDir =
                (outputDir != null && !outputDir.isEmpty())
                        ? outputDir
                        //            : LOCAL_WORKING_DIR + "toy.fasta";
                        //            : LOCAL_WORKING_DIR + "debug_fasta_subset.fasta";
                        : LOCAL_WORKING_DIR + "output/";
        JSONArray searchParameters =
                (serializedParameters != null && !serializedParameters.isEmpty())
                        ? new JSONArray(serializedParameters)
                        : new JSONArray(
                        "[{\"MS1_PPM\":\"30\",\"MS2_PPM\":\"10\",\"RT_TOL\":\"0.16\",\"MAX_PPT\":\"5\",\"MIN_PEPTIDE_LEN\":\"7\",\"PRECURSOR_Z\":\"1\",\"FRGION_Z\":\"2\"}]");

        // Sanity check for file existence(s)
        File mzXMLIn = new File(mzXMLFile);
        if (!mzXMLIn.exists()) {
            parser.printHelp();
            System.exit(1);
        }
        rawFileName = Paths.get(mzXMLFile).getFileName().toString();

        System.out.println("Configuring to extract MS2 XIC & FASTA Match ...");

        // MS2 XIC Extraction Outfile
        String xicOutFile = outDir + "isolationWindowRanges.out";

        // TODO: Have a more elegant way; of parsing the params
        JSONObject obj = searchParameters.getJSONObject(0);
        Map<String, Object> parameterMap = obj.toMap();
        Object ms2_ppm_object = parameterMap.get("MS2_PPM");
        int ms2_ppm;
        if (ms2_ppm_object instanceof String) {
            ms2_ppm = Integer.parseInt((String) ms2_ppm_object);
        } else {
            ms2_ppm = (int) parameterMap.get("MS2_PPM");
        }

        // Configure for timing the matching process
        StopWatch fullRunStopWatch = new StopWatch();
        fullRunStopWatch.start();

//    /* MS1 Feature Detection
        //* TODO: Uncomment for Xiangyuan's MS1 Feature Detection
        // Step 1: Identify possible precursors using Xiangyuan's MS1 component
        StopWatch ms1FeatureDetectionStopWatch = new StopWatch();
        ms1FeatureDetectionStopWatch.start();

        StopWatch ms1FeatureReadStopWatch = new StopWatch();
        ms1FeatureReadStopWatch.start();
        OpenMzxml of = new OpenMzxml(mzXMLFile);

        ms1FeatureReadStopWatch.stop();
        System.out.println("Read  mzML File Elapsed Time in Minutes: " + ms1FeatureReadStopWatch.getTime(TimeUnit.MINUTES) + " ...");


        FeatureDetect featureDetect = new FeatureDetect(of, FeatureDetect.DetectionType.MS1);
        featureDetect.detectFeatures(mzXMLFile);

        // write to mzXMLFile_svr_score
        System.out.println("Running SVR.py ...");
        ProcessBuilder processBuilder = new ProcessBuilder();
        ArrayList<String> svrCommand = new ArrayList<>(processBuilderCommand);
        svrCommand.add("/waterlooms/src/main/python/SVR.py");
        svrCommand.add("-features");
        svrCommand.add(DOCKER_WORKING_DIR + rawFileName);
        processBuilder.command(svrCommand);
        try {
            // TODO: Suppress output if necessary for python
            Process process = processBuilder.inheritIO().start();
            int exitCode = process.waitFor();
            System.out.println("\nSVR.py exited with code : " + exitCode);
        } catch (InterruptedException $e) {
            $e.printStackTrace();
        }

        // write to mzXMLFile_feature_one_z
        System.out.println("Selecting charge with FeatureSelect ...");
        FeatureSelect featureSelect = new FeatureSelect();
        featureSelect.selectFeature(
                LOCAL_WORKING_DIR + rawFileName);

        // Step 4 : write to mzXMLFile_nn_score
        ArrayList<String> nnCommand = new ArrayList<>(processBuilderCommand);
        nnCommand.add("/waterlooms/src/main/python/NN.py");
        nnCommand.add("-features");
        nnCommand.add(DOCKER_WORKING_DIR + rawFileName);

        processBuilder.command(nnCommand);
        try {
            Process process = processBuilder.start();
            int exitCode = process.waitFor();
            System.out.println("\nNN.py exited with code : " + exitCode);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

        // Step 5 : write to mzXMLFile_feature
        // Updated to just save as a list of SVRScores
        featureSelect.finalizeFeatureToMS1(
                LOCAL_WORKING_DIR + rawFileName,
                LOCAL_WORKING_DIR, rawFileName
        );

//    add isotope single with single trail
        featureDetect.searchIsotopeWithSignalToMS1(featureSelect.getSpec(), 4, 8);
        //combine and remove duplicate trail
        featureSelect.getSpec().combinRaw();

        System.out.println("Completed MSTracer Precursor Detection for " + rawFileName + "...");
        ms1FeatureDetectionStopWatch.stop();
        System.out.println("MS1 Precursor Detection Elapsed Time in Minutes: " + ms1FeatureDetectionStopWatch.getTime(TimeUnit.MINUTES) + " ...");

//    /* MS2 TrailDetection
        // Re-Parse the saved precursor list and order the precursors by RT
        //* STEP 2: Detect MS2 Trails
        StopWatch ms2FeatureDetectionStopWatch = new StopWatch();
        ms2FeatureDetectionStopWatch.start();

        TreeMap<Double, List<MSOneTrail>> mapRTMSoneTrail = featureSelect.getSpec().arrMSOneTrail.stream().sorted()
                .collect(Collectors.groupingBy(MSOneTrail::getRt, TreeMap::new, Collectors.toList()));

        GenerateMS2Trail genMS2 = new GenerateMS2Trail();

        genMS2.featureDetectShouldReturnAnArrayListOfXICs(filepath, of, mapRTMSoneTrail);//directory and mzXML filename

//
//    IsolationWindow.writeIsolationWindowPair(isolation_windows, xicOutFile);
//    System.out.println("Writing of XIC .ser and .tsv completed!");
        ms2FeatureDetectionStopWatch.stop();
        System.out.println("MS2 Precursor Detection Elapsed Time in Minutes: " + ms2FeatureDetectionStopWatch.getTime(TimeUnit.MINUTES) + " ...");
        // */


        // read fasta
        StopWatch peptideGenerateStopWatch = new StopWatch();
        peptideGenerateStopWatch.start();
        ArrayList<edu.uw.waterlooms.peptideMatch.Genome> genomes = edu.uw.waterlooms.peptideMatch.FastaFile.ReadFile(fastaFile);

        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top

        List<edu.uw.waterlooms.peptideMatch.Peptide> pepsOri = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= edu.uw.waterlooms.peptideMatch.Utils.thresholdPeptideSizeMin).sorted().collect(Collectors.toCollection(ArrayList::new));

        if (edu.uw.waterlooms.peptideMatch.Utils.bRemovePeptideFromBothTargetAndDecoy) {

            TreeMap<String, List<edu.uw.waterlooms.peptideMatch.Peptide>> mapPepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(edu.uw.waterlooms.peptideMatch.Peptide::getComposition, TreeMap::new, Collectors.toList()));

            TreeMap<String, List<edu.uw.waterlooms.peptideMatch.Peptide>> mapILreplacePepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(edu.uw.waterlooms.peptideMatch.Peptide::getReplaceILComposition, TreeMap::new, Collectors.toList()));

            System.out.println(mapPepstrListPep.size());
            ArrayList<edu.uw.waterlooms.peptideMatch.Peptide> arrPepsRemove = new ArrayList<>();

            for (String strPep : mapPepstrListPep.keySet()) {
                boolean bDecoyType = false;
                boolean bTargetType = false;

                if (mapPepstrListPep.get(strPep).size() > 1) {
                    for (edu.uw.waterlooms.peptideMatch.Peptide pep : mapPepstrListPep.get(strPep)) {
                        boolean isDecoy = pep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }
                }
                //REMOVE I L MUTATE peptide in both decoy and target
                if (strPep.contains("I") || strPep.contains("L")) {
                    for (edu.uw.waterlooms.peptideMatch.Peptide muSamePep : mapILreplacePepstrListPep.get(strPep.replaceAll("L", "I"))) {
                        boolean isDecoy = muSamePep.id.contains("DeBruijn");
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
        List<edu.uw.waterlooms.peptideMatch.Peptide> peps = pepsOri.stream().distinct().sorted().collect(Collectors.toCollection(ArrayList::new));

        double[] pepsMass = peps.stream().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());
        peptideGenerateStopWatch.stop();
        System.out.println("Peptide Generate Elapsed Time in Minutes: " + peptideGenerateStopWatch.getTime(TimeUnit.MINUTES) + " ...");


        String psmOutfile = (LOCAL_WORKING_DIR + rawFileName).replaceFirst("[.][^.]+$", "");

        FileWriter fileWriter = new FileWriter(psmOutfile + "_psm.csv");
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tmsoneMZ\tmsoneZ\tmsonePeakAreaLocalRank" +
                "\tms2rt\tpeptideLength\tMS1PeptideMassError\tdPeptideMutationRate\tibyComplemetaryCount\tdCosWithPredictMSMS\tdMS1WithMS2TrailCoscinSimilarity\tMS2TrailCoscinSimilarity\tBIon\tYIon");
        for (int iw = 0; iw < Utils.iMS2Charge; iw++) {
            bufferedWriter.write("\tdlogBionPeakAreaSum_C" + (iw + 1) + "\tdlogYionPeakAreaSum_C" + (iw + 1) + "\tdBionMassErrorSum_C" + (iw + 1) + "\tdYionMassErrorSum_C" + (iw + 1) + "\tiBionPeakHalfSum_C" + (iw + 1) + "\tiYionPeakHalfSum_C" + (iw + 1) + "\tbionMatchCount_C" + (iw + 1) + "\tyionMatchCount" +
                    "_C" + (iw + 1) + "\tdBionRetentionTimeErrorSum_C" + (iw + 1) + "\tdYionRetentionTimeErrorSum_C" + (iw + 1) + "\tiBConsective_C" + (iw + 1) + "\tiYConsective_C" + (iw + 1) + "\tdBionCosWithPredictMSMS_C" + (iw + 1) + "\tdYionCosWithPredictMSMS_C" + (iw + 1));
        }
        bufferedWriter.write("\tarrMatchedBion\tarrMatchedYion\tadjustBYIon\tSoredIndex\twindowSize\tdBionCos60WithPredictMSMS_C1\tdCosWithPredictMSMSMatchedPredict\tprositRT" +
                "\tibAllMatchedWithloss\tibMatchNormalAndH2Oloss\tibMatchNormalAndNH3loss\tiyAllMatchedWithloss\tiyMatchNormalAndH2Oloss\tiyMatchNormalAndNH3loss\tisRawdata" +
                "\tpeptideAnomiAcidFreRate\tlastIonNumBC1\tlastIonNumYC1\ttop6allCos\ttop6c1Cos\ttop6Bc1Cos\ttop6Yc1Cos\tcount"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
                + "\n");


        //output multiple file, concise file, only peptidefile
        FileWriter concisefileWriter;
        BufferedWriter concisebufferedWriter = null;
        FileWriter peptidesForpDeep3Writer = null;
        BufferedWriter peptidesForpDeep3bufferedWriter = null;
        FileWriter peptidesForpPrositWriter = null;
        BufferedWriter peptidesForpPrositbufferedWriter = null;
        if (edu.uw.waterlooms.peptideMatch.Utils.bzOutputConciseAndPeptideFiles) {
            concisefileWriter = new FileWriter( LOCAL_WORKING_DIR +"concisefile_0502.csv");
            concisebufferedWriter = new BufferedWriter(concisefileWriter);
            peptidesForpDeep3Writer = new FileWriter( LOCAL_WORKING_DIR +"peptidesForpDeep3_0502.csv");
            peptidesForpDeep3bufferedWriter = new BufferedWriter(peptidesForpDeep3Writer);
            peptidesForpPrositWriter = new FileWriter( LOCAL_WORKING_DIR +"peptidesForpProsit_0502.csv");
            peptidesForpPrositbufferedWriter = new BufferedWriter(peptidesForpPrositWriter);
            concisebufferedWriter.write("msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tgetQuality_score\tmsonemz\tmsonecharge\tBIon\tYIon" + "\n");
            peptidesForpDeep3bufferedWriter.write("peptide\tmodinfo\tcharge" + "\n");
            peptidesForpPrositbufferedWriter.write("modified_sequence\tcollision_energy\tprecursor_charge\tfragmentation\n");

        }
//get ms2 intensity and rt predicted by prosit
        StopWatch prositResultStopWatch = new StopWatch();
        prositResultStopWatch.start();

//        Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit(LOCAL_WORKING_DIR);//firstrunn 0829
        Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit(prositDir);//firstrunn 1123 //要求输入prositResult 目录

        prositResultStopWatch.stop();
        System.out.println("MS2 Intensity and RT predicted by prosit Elapsed Time in Minutes: " + prositResultStopWatch.getTime(TimeUnit.MINUTES) + " ...");

        System.out.println("Matching Peptides to MS2 XICs ...");


        featureSelect.getSpec().clearRTSandMapRT();


        double[] arrRts = featureSelect.getSpec().arrMSOneTrail.stream().mapToDouble(MSOneTrail::getRt)
                .distinct()
                .sorted()
                .toArray();

        edu.uw.waterlooms.peptideMatch.IsolationWindowCollection iswc = genMS2.iswc;

        //set ms1feature peakarea rank

        AtomicInteger ims1 = new AtomicInteger();


        ProgressBar pbMS1 = new ProgressBar("Progress MS1 area rank", arrRts.length);
        pbMS1.start();


        featureSelect.getSpec().arrMSOneTrail.stream().parallel().forEach(ms1 ->
        {
            MSOneTrail ms1Trail = ms1;
            pbMS1.step();
            long iPeakAreaRank = 1;//至少排名第一
            List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(ms1Trail.getMz(), 1, true);
            for (MSTwoTrailSet window : arrwindow) {
                //get ms1feature RT distance between Math.abs(y-v) <= Utils.thresholdRT
                int iMS1RTPosBegin = Arrays.binarySearch(arrRts, ms1Trail.getRt() - edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1CheckPeakAreaRT);//正负15秒内
                iMS1RTPosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);

                for (int iRTPos = iMS1RTPosBegin;
                     iRTPos < arrRts.length && arrRts[iRTPos] < ms1Trail.getRt() + edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1CheckPeakAreaRT;
                     iRTPos++) {
                    if (window != null) {
                        //0914 multi thread
                        iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
                                x.getZ() == ms1Trail.getZ() && x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh &&
                                        x.getQuantification_peaks_area() > ms1Trail.getQuantification_peaks_area()
                        ).count();
                    }
                }
            }
            ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));

        });
        pbMS1.stop();

        Map<Double,Set<String>> mapRtPepSet = Collections.synchronizedMap(new HashMap<>());
        for(double drt:arrRts)
        {
            Set<String> pepSet = Collections.synchronizedSet(new HashSet<>());
            mapRtPepSet.put(drt,pepSet);
        }
        List<MatchPeptideWithMultiCharge> psmTop1_Top6w = new LinkedList<>();

        //为节省内存，一个一个Window process
//        for(int iw=0;iw<iswc.windows.size();iw++)
//        int i=0;
        for(MSTwoTrailSet ms2window:iswc.windows)
        {
//test
//            i++;
//            if(i>2) {
//                break;
//
//            }
            List<MatchPeptideWithMultiCharge> listPeptideAdd     = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());
            List<MatchPeptideWithMultiCharge> listPeptideAddTrue = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());
            AtomicInteger icount = new AtomicInteger();


            ProgressBar pb = new ProgressBar("Progress", arrRts.length);
            pb.start();
            Arrays.stream(arrRts).parallel().forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();
                        pb.step();
                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);

//                        Set<String> uniquePeptideString = Collections.synchronizedSet(new HashSet<>());
//                        Set<String> uniquePeptideString = mapRtPepSet.get(y);


                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {
                                        //是否在当前窗口
                                        if (msonetrail.getMz()<=ms2window.mzHigh && msonetrail.getMz()>=ms2window.mzLow &&
                                                msonetrail.getMass() - edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {
                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window

                                            LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide = new LinkedList<MatchPeptideWithMultiCharge>();
                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {
                                                //begin of window with mz of ms1 feature
                                                // if mass matches, then pep match
//                                                List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(msonetrail.getMz(), 1, true);
//                                                for (MSTwoTrailSet window : arrwindow) {
//                                                    MS1MatchPeptide_top6PredictIntensity(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msonetrail, listMaxiumPeptide, peps.get(iPos), window);
                                                    MS1MatchPeptide_top6PredictIntensity(icount, peps, bufferedWriter,  MapMS2Intensity, msonetrail, listMaxiumPeptide, peps.get(iPos), ms2window);
//                                                }
                                            }

//                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);

                                            if (listMaxiumPeptide.size() > 0) {
                                                if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputALL) {
                                                    listPeptideAddTrue.addAll(listMaxiumPeptide);

                                                } else if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputMax) {
                                                    listPeptideAddTrue.add(listMaxiumPeptide.get(0));
                                                } else if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputSecondMax) {
                                                    if (listMaxiumPeptide.size() > 1) {
                                                        listPeptideAddTrue.add(listMaxiumPeptide.get(1));
                                                    }
                                                }
                                            }
//                                            listPeptideAddTrue.addAll(listMaxiumPeptide);
                                            ims1.addAndGet(listMaxiumPeptide.size());
                                            if (edu.uw.waterlooms.peptideMatch.Utils.bzCheckPeptideWithAllChargeAtSameRT) {
                                                //检查该MS1相同RT且已经匹配的相同PEPTIDE里面所有电荷的二级谱
                                                int iCurrenMS1Z = msonetrail.getZ();
//                                                Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());
                                                for (MatchPeptideWithMultiCharge matchPep : listMaxiumPeptide) {
                                                    //debug
//                                                    if(matchPep.pep.composition.equals("GGDRVPADIR"))
//                                                    {
//                                                        System.out.println("ms1rt:"+y+"      "+JSONWriter.valueToString(mapRtPepSet.get(y)));
//                                                    }
                                                    //是否重复出现的PEPTIDE，若相同只进来一次
//                                                    if (uniquePeptideString.contains(matchPep.pep.composition)) {
                                                    if (mapRtPepSet.get(y).contains(matchPep.pep.composition)) {
                                                        continue;
                                                    }
                                                    for (int iAllC = 1; iAllC <= edu.uw.waterlooms.peptideMatch.Utils.iMS1Charget; iAllC++) {
                                                        //condition is the z is larger than current z
                                                        if (iAllC == iCurrenMS1Z) {
                                                            continue;//means this
                                                        }
                                                        double dPepMZWithC = edu.uw.waterlooms.peptideMatch.Utils.MassToMz(matchPep.pep.mass, iAllC);
                                                        //search the mz in lMsoneTrail,如果存在（可能存在一个或多个）且比当前小就不search 二级谱，否则search，然后在标记该peptide已经被search，下次就不再search

                                                        if (iAllC < iCurrenMS1Z) {
                                                            //if it is matched, it should told the problem there is no matched at successive process(create a map??)
                                                            //get the window
                                                            boolean isE = false;
                                                            int iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC * (1 - edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            for (int iPos = iMS1PosBegin;
                                                                 iPos < lMsoneTrail.size() && lMsoneTrail.get(iPos).getMz() < dPepMZWithC * (1 + edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM);//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                                 iPos++) {

                                                                isE = true;
                                                                break;

                                                            }
                                                            if (isE) {
                                                                continue;
                                                            }

                                                        } else if (iAllC > iCurrenMS1Z) {
                                                            //(the larger z means mz is less than current mz that has been checked before)
                                                            int iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC * (1 - edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            if (iMS1PosBegin
                                                                    < lMsoneTrail.size() && lMsoneTrail.get(iMS1PosBegin).getMz() < dPepMZWithC * (1 + edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM))//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                            {

                                                                continue;
                                                            }


                                                        }
//                                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(dPepMZWithC, 1);
                                                        List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(dPepMZWithC, 1, true);
                                                        for (MSTwoTrailSet window : arrwindow) {
                                                            MSOneTrail msOneTrail = new MSOneTrail(-1, dPepMZWithC, msonetrail.getRt(), iAllC, 0, 0, 0);
                                                            LinkedList<MatchPeptideWithMultiCharge> listPeptide = new LinkedList<MatchPeptideWithMultiCharge>();
//                                                            MS1MatchPeptide(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msOneTrail, listPeptide, matchPep.pep, window);
                                                            MS1MatchPeptide(icount, peps, bufferedWriter,  MapMS2Intensity, msOneTrail, listPeptide, matchPep.pep, window);
//                                                            if(matchPep.pep.composition.equals("GGDRVPADIR"))
//                                                            {
//                                                                System.out.println("ms1rt:"+y+"     GGDRVPADIR ");
//                                                            }
//
//                                                            可以做一些条件限制，节省空间，并不是每一个虚拟的precursor所对应的prsm都用到
//                                                            msonefeature == -1 &
//                                                                    matchscore > 1 & (MS2TrailCoscinSimilarity / peptideLength) > 0.40)
                                                            for(MatchPeptideWithMultiCharge matchPSM:listPeptide)
                                                            {
                                                                if(matchPSM.dMatch>Utils.dThresholdVirturalPrecursorPSMMatchScore)
//                                                                        && matchPSM.pep.composition.length()>Utils.dThresholdVirturalPrecursorPSMCoscinSimilarityScore)
                                                                {
                                                                    listPeptideAdd.add(matchPSM);


                                                                }
                                                            }
//                                                            listPeptideAdd.addAll(listPeptide);
                                                        }

                                                        //check other charge
                                                    }
                                                    //check whether it has matched before

//                                                    uniquePeptideString.add(matchPep.pep.composition);
                                                    mapRtPepSet.get(y).add(matchPep.pep.composition);
//                                                    if(matchPep.pep.composition.equals("GGDRVPADIR"))
//                                                    {
//                                                        System.out.println("ms1rt:"+y+"      "+JSONWriter.valueToString(mapRtPepSet.get(y)));
//                                                    }
                                                }


                                            }
                                        }

                                    }
                            );
                        }
                    }
            );
            pb.stop();


            System.out.println("add listPeptideAddTrue size:" + listPeptideAddTrue.size());

            System.out.println("add CheckPeptideWithAllChargeAtSameRT size:" + listPeptideAdd.size());
            List<MatchPeptideWithMultiCharge> sorted = listPeptideAdd;
            sorted.addAll(listPeptideAddTrue);

            System.out.println("Original Count :" + ims1);

            System.out.println("after combine CheckPeptideWithAllChargeAtSameRT size:" + sorted.size());


            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long, Integer> mapIDCount = new HashMap();


            Collections.sort(sorted, Collections.reverseOrder());


            System.out.println("---------OutputALL Begin---------");
            ProgressBar pboutput = new ProgressBar("OutputALL", sorted.size());
            pboutput.start();
            int iTop6w = 0;
            for (MatchPeptideWithMultiCharge mp : sorted) {
                pboutput.step();
                if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        if(iTop6w<Utils.iCount1TopPSM)
                        {
                            psmTop1_Top6w.add(mp);
                        }
                        iTop6w++;
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputALL) {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");

                } else {
                    bufferedWriter.write(mp.toString() + "\n");


                }
                if (edu.uw.waterlooms.peptideMatch.Utils.bzOutputConciseAndPeptideFiles) {

                    concisebufferedWriter.write(mp.getConciseString() + "\n");
                    peptidesForpDeep3bufferedWriter.write(mp.pep.composition + "\t\t" + mp.msOneTrail.getZ() + "\n");
                    peptidesForpPrositbufferedWriter.write(mp.pep.composition + "\t27.5\t" + mp.msOneTrail.getZ() + "\tHCD\n");

                }

                iCountForFDR++;

            }

            //keep top N elements
            psmTop1_Top6w =  psmTop1_Top6w.stream().sorted(Collections.reverseOrder())
                    .limit(Utils.iCount1TopPSM)
                    .collect(Collectors.toList());


            //for gc
            sorted=null;
            mapIDCount=null;

            pboutput.stop();
            System.gc();

        }

        System.out.println("---------OutputALL Finish---------");

        if (Utils.bzOutputConciseAndPeptideFiles) {
            concisebufferedWriter.flush();
            peptidesForpDeep3bufferedWriter.flush();
            peptidesForpPrositbufferedWriter.flush();
            concisebufferedWriter.close();
            peptidesForpDeep3bufferedWriter.close();
            peptidesForpPrositbufferedWriter.close();
        }


        bufferedWriter.flush();
        bufferedWriter.close();


//for autort rt and irt alignment

        String autoRtOut = autortFile;//"autoRTMouseR021t3d.tsv";
        Map<String,Double> mapAutoRT = new HashMap<>();
        FileReader freaderAutort;
        try {
            freaderAutort = new FileReader(DOCKER_WORKING_DIR + autoRtOut);
            try (BufferedReader br = new BufferedReader(freaderAutort)) {
                String line;
                br.readLine();//skip first line
                while ((line = br.readLine()) != null) {
                    String[] range = line.split("\\s+");
                    mapAutoRT.put(range[0],Double.parseDouble(range[2]));
                }
            }

        } catch (IOException noFile) {
            throw new FileNotFoundException();
        }
        FileWriter fileWriterTop6w = new FileWriter(psmOutfile + "_psmTop6w.csv");
        BufferedWriter bufferedWriterTop6w = new BufferedWriter(fileWriterTop6w);
        bufferedWriterTop6w.write("RT\tiRT\n");
        for(MatchPeptideWithMultiCharge mpT:psmTop1_Top6w)
            bufferedWriterTop6w.write(mpT.msOneTrail.getRt() +"\t" + mapAutoRT.get(mpT.pep.composition)+ "\n");

        bufferedWriterTop6w.flush();
        bufferedWriterTop6w.close();
        fileWriterTop6w.close();



//20221117
        // write to mzXMLFile_svr_score
        System.out.println("Running rt.py ...");
        ProcessBuilder processBuilderRT = new ProcessBuilder();
        ArrayList<String> rtCommand = new ArrayList<>(processBuilderCommand);
        rtCommand.add("/waterlooms/src/main/python/rt_utils.py");
        rtCommand.add("-features");
//        rtCommand.add(DOCKER_WORKING_DIR + "humanR01forRTminfit_maxfit_diff_1t3d.csv");
        rtCommand.add(psmOutfile+ "_psmTop6w.csv");

        System.out.println(rtCommand);
        processBuilderRT.command(rtCommand);
        try {
            // TODO: Suppress output if necessary for python
            Process process = processBuilderRT.inheritIO().start();
            int exitCode = process.waitFor();
            System.out.println("\nRT.py exited with code : " + exitCode);
        } catch (InterruptedException $e) {
            $e.printStackTrace();
        }



        // test another docker and run r code
//        ArrayList<String> processBuilderRCommand = new ArrayList<String>();
//        processBuilderRCommand.add("docker");
//        processBuilderRCommand.add("exec");
//        processBuilderRCommand.add("-t");
//        processBuilderRCommand.add("rmodel_dp");//need change
//
//        processBuilderRCommand.add("/usr/local/lib/R/bin/Rscript");
//
//        System.out.println("Running model.r ...");
//        ProcessBuilder processBuilderRL = new ProcessBuilder();
//        ArrayList<String> rlCommand = new ArrayList<>(processBuilderRCommand);
//        rlCommand.add("/home/rstudio/data/model.R" );


        // test another docker and run r code
        ArrayList<String> processBuilderRCommand = new ArrayList<String>();
        processBuilderRCommand.add("Rscript");

        System.out.println("Running model.r ...");
        ProcessBuilder processBuilderRL = new ProcessBuilder();
        ArrayList<String> rlCommand = new ArrayList<>(processBuilderRCommand);
        rlCommand.add("/waterlooms/src/main/R/model.R" );
        rlCommand.add(rawFileName.replaceFirst("[.][^.]+$", ""));//FileName
        rlCommand.add(psmOutfile+ "_psm.csv");//database searching result
        rlCommand.add(autortFile);//("autoRTHumanR01_1t3d.csv");//autoRTresult
        rlCommand.add(psmOutfile + "_psmTop6w_percentile_list.csv");//("humanR01_RTminfit_maxfit_diff_1t3d.txt");//RT alignment
        rlCommand.add(DOCKER_WORKING_DIR);
        processBuilderRL.command(rlCommand);
        try {
            // TODO: Suppress output if necessary for python
            Process process = processBuilderRL.inheritIO().start();
            int exitCode = process.waitFor();
            System.out.println("\nmodel.r exited with code : " + exitCode);
        } catch (InterruptedException $e) {
            $e.printStackTrace();
        }

//        docker exec -t lucid_hawking /usr/local/lib/R/bin/Rscript /home/rstudio/data/model.R

//        System.exit(1);

//        rlCommand.add(rawFileName.replaceFirst("[.][^.]+$", ""));
//        rlCommand.add(psmOutfile+ "_psm.csv");
//        rlCommand.add(autoRtOut);
//        rlCommand.add(psmOutfile+ "_percentile_list.csv");
//        processBuilderRL.command(rlCommand);
//        try {
//            // TODO: Suppress output if necessary for python
//            Process process = processBuilderRL.inheritIO().start();
//            int exitCode = process.waitFor();
//            System.out.println("\nmodel.r exited with code : " + exitCode);
//        } catch (InterruptedException $e) {
//            $e.printStackTrace();
//        }

//        // test another docker and run r code
//        ArrayList<String> processBuilderRCommand = new ArrayList<String>();
//        processBuilderRCommand.add("docker");
//        processBuilderRCommand.add("exec");
//        processBuilderRCommand.add("-t");
//        processBuilderRCommand.add("lucid_hawking");//need change
//
//        processBuilderRCommand.add("/usr/local/lib/R/bin/Rscript");
//
//        System.out.println("Running model.r ...");
//        ProcessBuilder processBuilderRL = new ProcessBuilder();
//        ArrayList<String> rlCommand = new ArrayList<>(processBuilderRCommand);
//        rlCommand.add("/home/rstudio/data/model.R" );
//        rlCommand.add(psmOutfile+"_prsm.csv");
//        rlCommand.add(DOCKER_WORKING_DIR+"autort.csv");
//        rlCommand.add(psmOutfile + "_psmTop6w_percentile_list.csv");
//        processBuilderRL.command(rlCommand);
//        try {
//            // TODO: Suppress output if necessary for python
//            Process process = processBuilderRL.inheritIO().start();
//            int exitCode = process.waitFor();
//            System.out.println("\nmodel.r exited with code : " + exitCode);
//        } catch (InterruptedException $e) {
//            $e.printStackTrace();
//        }

//        docker exec -t lucid_hawking /usr/local/lib/R/bin/Rscript /home/rstudio/data/model.R

//        System.exit(1);

        System.out.println("Matched with Candidate Peptides ...");

        fullRunStopWatch.stop();
        System.out.println("Elapsed Time in Minutes: " + fullRunStopWatch.getTime(TimeUnit.MINUTES) + " ...");
        System.exit(0);
    }

    public static void main_readfromFile(String[] args) throws IOException {
        // Configure for timing the matching process
        StopWatch fullRunStopWatch = new StopWatch();
        fullRunStopWatch.start();
        // Parse environment variable set in Dockerfile to validate if the executable is running within
        // the container
        boolean containerExecution = Boolean.parseBoolean(System.getenv("RUN_WITHIN_CONTAINER"));
        ArrayList<String> processBuilderCommand = new ArrayList<String>();
        if (!containerExecution) {
            processBuilderCommand.add("docker-compose");
            processBuilderCommand.add("exec");
            processBuilderCommand.add("-T");
            processBuilderCommand.add("waterlooms");
        }
        processBuilderCommand.add("python");

        Path path = FileSystems.getDefault().getPath("").toAbsolutePath();
        String DOCKER_WORKING_DIR = "/data/";
        String LOCAL_WORKING_DIR = containerExecution ? "/data/" : path + "/data/";

        ArgumentParser parser =
                ArgumentParsers.newFor("java -jar waterlooms.jar")
                        .build()
                        .defaultHelp(true)
                        .description("DIA Analysis.");
        parser
                .addArgument("-mzXML")
                .metavar("FILE")
                .type(Arguments.fileType().acceptSystemIn().verifyExists())
                .help("mzXML File");
        parser.addArgument("-parameters").type(String.class).help("JSON String Dict of Parameters");
        parser
                .addArgument("-fasta")
                .metavar("FILE")
                .type(Arguments.fileType().acceptSystemIn().verifyExists())
                .help("Fasta Peptide Database File. Must contain decoys.");
        parser
                .addArgument("-outputDir")
                .type(String.class)
                .help("Directory for outputting the result CSV.");

        // Parse Arguments
        Namespace ns = null;
        try {
            ns = parser.parseArgs(args);
        } catch (ArgumentParserException $exception) {
            parser.handleError($exception);
            System.exit(1);
        }

        String mzXMLInFile = ns.getString("mzXML");
        String fastaInFile = ns.getString("fasta");
        String outputDir = ns.getString("outputDir");
        String serializedParameters = ns.getString("parameters");

        // Configure the file for local debugging
//    String rawFileName = "fig2_mp_dia_multi_organism_run_mhrm_r01.mzXML";
//    String rawFileName = "fig2_mp_dia_multi_organism_run_mhrm_r01.mzXML";
        //   String rawFileName = "r01_dia_data_deconvoluted_noms1.mzXML";
        String rawFileName = "toy.mzXML";


        if ((mzXMLInFile != null && !mzXMLInFile.isEmpty())) {
            // TODO: DIA-WEBAPP submits mzXMLInFile as a path/file.mzXML
            // TODO: Need to strip this and set rawFileName
            rawFileName = mzXMLInFile;
        }


        // Set defaults from main/resources/ folder if arguments not specified (development purposes)
        String mzXMLFile =
                (mzXMLInFile != null && !mzXMLInFile.isEmpty())
                        ? mzXMLInFile
                        : LOCAL_WORKING_DIR + rawFileName;


        //follow mzxmlFile path set the docker_working_path
        String filepath = mzXMLFile.substring(0, mzXMLFile.lastIndexOf("/"));
        DOCKER_WORKING_DIR = filepath + "/";
        LOCAL_WORKING_DIR = filepath + "/";
        System.out.println(DOCKER_WORKING_DIR);

        String fastaFile =
                (fastaInFile != null && !fastaInFile.isEmpty())
                        ? fastaInFile
//            : LOCAL_WORKING_DIR + "database_withdecoy.fasta";
                        : LOCAL_WORKING_DIR + "homo_sapiens_withdecoy.fasta";
        String outDir =
                (outputDir != null && !outputDir.isEmpty())
                        ? outputDir
                        //            : LOCAL_WORKING_DIR + "toy.fasta";
                        //            : LOCAL_WORKING_DIR + "debug_fasta_subset.fasta";
                        : LOCAL_WORKING_DIR + "output/";



        MSOneTrailSet spec = new MSOneTrailSet();
        if(!Utils.bzOnlyRawdata)
        {
            spec.readFromeTrailFile(Config.spectrumFolder +Utils.strinputPrecursor,false);//读取xiangyuan的MS1 FEATURE文件
            spec.generateRTSandMapRT();

        }
        spec.readFromeTrailFile(Config.spectrumFolder + Utils.strRawIsotopeFileName,true);//读取raw的isotope文件
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
        edu.uw.waterlooms.peptideMatch.IsolationWindowCollection iswc = new edu.uw.waterlooms.peptideMatch.IsolationWindowCollection();
        iswc.IsolationWindowCollection_paralle(Config.spectrumFolder + Utils.strIsolationWindowCollectionFileName,mapRTMSoneTrail);



        // read fasta
        StopWatch peptideGenerateStopWatch = new StopWatch();
        peptideGenerateStopWatch.start();
        ArrayList<edu.uw.waterlooms.peptideMatch.Genome> genomes = edu.uw.waterlooms.peptideMatch.FastaFile.ReadFile(fastaFile);

        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top

        List<edu.uw.waterlooms.peptideMatch.Peptide> pepsOri = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= edu.uw.waterlooms.peptideMatch.Utils.thresholdPeptideSizeMin).sorted().collect(Collectors.toCollection(ArrayList::new));

        if (edu.uw.waterlooms.peptideMatch.Utils.bRemovePeptideFromBothTargetAndDecoy) {

            TreeMap<String, List<edu.uw.waterlooms.peptideMatch.Peptide>> mapPepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(edu.uw.waterlooms.peptideMatch.Peptide::getComposition, TreeMap::new, Collectors.toList()));

            TreeMap<String, List<edu.uw.waterlooms.peptideMatch.Peptide>> mapILreplacePepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(edu.uw.waterlooms.peptideMatch.Peptide::getReplaceILComposition, TreeMap::new, Collectors.toList()));

            System.out.println(mapPepstrListPep.size());
            ArrayList<edu.uw.waterlooms.peptideMatch.Peptide> arrPepsRemove = new ArrayList<>();

            for (String strPep : mapPepstrListPep.keySet()) {
                boolean bDecoyType = false;
                boolean bTargetType = false;

                if (mapPepstrListPep.get(strPep).size() > 1) {
                    for (edu.uw.waterlooms.peptideMatch.Peptide pep : mapPepstrListPep.get(strPep)) {
                        boolean isDecoy = pep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }
                }
                //REMOVE I L MUTATE peptide in both decoy and target
                if (strPep.contains("I") || strPep.contains("L")) {
                    for (edu.uw.waterlooms.peptideMatch.Peptide muSamePep : mapILreplacePepstrListPep.get(strPep.replaceAll("L", "I"))) {
                        boolean isDecoy = muSamePep.id.contains("DeBruijn");
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
        List<edu.uw.waterlooms.peptideMatch.Peptide> peps = pepsOri.stream().distinct().sorted().collect(Collectors.toCollection(ArrayList::new));

        double[] pepsMass = peps.stream().mapToDouble(x -> x.mass).toArray();

        System.out.println("The number of peps is: " + peps.size());
        peptideGenerateStopWatch.stop();
        System.out.println("Peptide Generate Elapsed Time in Minutes: " + peptideGenerateStopWatch.getTime(TimeUnit.MINUTES) + " ...");


        String psmOutfile = (LOCAL_WORKING_DIR + rawFileName).replaceFirst("[.][^.]+$", "");

        FileWriter fileWriter = new FileWriter(psmOutfile + "_psm.csv");
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

//        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmstwotime\tmatch\tmatchscore\tBIon\tYIon"+"\n");
        bufferedWriter.write("FDR\tmsonefeature\tmatchscore\tproteinname\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tmsoneQualityScore\tmsoneMZ\tmsoneZ\tmsonePeakAreaLocalRank" +
                "\tms2rt\tpeptideLength\tMS1PeptideMassError\tdPeptideMutationRate\tibyComplemetaryCount\tdCosWithPredictMSMS\tdMS1WithMS2TrailCoscinSimilarity\tMS2TrailCoscinSimilarity\tBIon\tYIon");
        for (int iw = 0; iw < Utils.iMS2Charge; iw++) {
            bufferedWriter.write("\tdlogBionPeakAreaSum_C" + (iw + 1) + "\tdlogYionPeakAreaSum_C" + (iw + 1) + "\tdBionMassErrorSum_C" + (iw + 1) + "\tdYionMassErrorSum_C" + (iw + 1) + "\tiBionPeakHalfSum_C" + (iw + 1) + "\tiYionPeakHalfSum_C" + (iw + 1) + "\tbionMatchCount_C" + (iw + 1) + "\tyionMatchCount" +
                    "_C" + (iw + 1) + "\tdBionRetentionTimeErrorSum_C" + (iw + 1) + "\tdYionRetentionTimeErrorSum_C" + (iw + 1) + "\tiBConsective_C" + (iw + 1) + "\tiYConsective_C" + (iw + 1) + "\tdBionCosWithPredictMSMS_C" + (iw + 1) + "\tdYionCosWithPredictMSMS_C" + (iw + 1));
        }
        bufferedWriter.write("\tarrMatchedBion\tarrMatchedYion\tadjustBYIon\tSoredIndex\twindowSize\tdBionCos60WithPredictMSMS_C1\tdCosWithPredictMSMSMatchedPredict\tprositRT" +
                "\tibAllMatchedWithloss\tibMatchNormalAndH2Oloss\tibMatchNormalAndNH3loss\tiyAllMatchedWithloss\tiyMatchNormalAndH2Oloss\tiyMatchNormalAndNH3loss\tisRawdata" +
                "\tpeptideAnomiAcidFreRate\tlastIonNumBC1\tlastIonNumYC1\ttop6allCos\ttop6c1Cos\ttop6Bc1Cos\ttop6Yc1Cos\tcount"//\tdBionCos60WithPredictMSMS_C2\tdYionCos60WithPredictMSMS_C1\tdYionCos60WithPredictMSMS_C2"/*	iBRepeatCount	iYRepeatCount	iBCoverCount	iYCoverCount	bStrCovered	yStrCovered*/
                + "\n");


        //output multiple file, concise file, only peptidefile
        FileWriter concisefileWriter;
        BufferedWriter concisebufferedWriter = null;
        FileWriter peptidesForpDeep3Writer = null;
        BufferedWriter peptidesForpDeep3bufferedWriter = null;
        FileWriter peptidesForpPrositWriter = null;
        BufferedWriter peptidesForpPrositbufferedWriter = null;
        if (edu.uw.waterlooms.peptideMatch.Utils.bzOutputConciseAndPeptideFiles) {
            concisefileWriter = new FileWriter( LOCAL_WORKING_DIR +"concisefile_0502.csv");
            concisebufferedWriter = new BufferedWriter(concisefileWriter);
            peptidesForpDeep3Writer = new FileWriter( LOCAL_WORKING_DIR +"peptidesForpDeep3_0502.csv");
            peptidesForpDeep3bufferedWriter = new BufferedWriter(peptidesForpDeep3Writer);
            peptidesForpPrositWriter = new FileWriter( LOCAL_WORKING_DIR +"peptidesForpProsit_0502.csv");
            peptidesForpPrositbufferedWriter = new BufferedWriter(peptidesForpPrositWriter);
            concisebufferedWriter.write("msonefeature\tdecoyTarget\tpeptide\tpepmass\tmsonemass\tmsonetime\tgetQuality_score\tmsonemz\tmsonecharge\tBIon\tYIon" + "\n");
            peptidesForpDeep3bufferedWriter.write("peptide\tmodinfo\tcharge" + "\n");
            peptidesForpPrositbufferedWriter.write("modified_sequence\tcollision_energy\tprecursor_charge\tfragmentation\n");

        }
//get ms2 intensity and rt predicted by prosit
        StopWatch prositResultStopWatch = new StopWatch();
        prositResultStopWatch.start();

        Map<String, MS2Intensity> MapMS2Intensity = MS2Intensity.getMapMS2IntensityFromProsit(LOCAL_WORKING_DIR);//firstrunn 0829

        prositResultStopWatch.stop();
        System.out.println("MS2 Intensity and RT predicted by prosit Elapsed Time in Minutes: " + prositResultStopWatch.getTime(TimeUnit.MINUTES) + " ...");

        System.out.println("Matching Peptides to MS2 XICs ...");


        spec.clearRTSandMapRT();





        //set ms1feature peakarea rank

        AtomicInteger ims1 = new AtomicInteger();


        ProgressBar pbMS1 = new ProgressBar("Progress MS1 area rank", arrRts.length);
        pbMS1.start();


        spec.arrMSOneTrail.stream().parallel().forEach(ms1 ->
        {
            MSOneTrail ms1Trail = ms1;
            pbMS1.step();
            long iPeakAreaRank = 1;//至少排名第一
            List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(ms1Trail.getMz(), 1, true);
            for (MSTwoTrailSet window : arrwindow) {
                //get ms1feature RT distance between Math.abs(y-v) <= Utils.thresholdRT
                int iMS1RTPosBegin = Arrays.binarySearch(arrRts, ms1Trail.getRt() - edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1CheckPeakAreaRT);//正负15秒内
                iMS1RTPosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(arrRts.length, iMS1RTPosBegin);

                for (int iRTPos = iMS1RTPosBegin;
                     iRTPos < arrRts.length && arrRts[iRTPos] < ms1Trail.getRt() + edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1CheckPeakAreaRT;
                     iRTPos++) {
                    if (window != null) {
                        //0914 multi thread
                        iPeakAreaRank += mapRTMSoneTrail.get(arrRts[iRTPos]).stream().filter((x) ->
                                x.getZ() == ms1Trail.getZ() && x.getMz() >= window.mzLow && x.getMz() <= window.mzHigh &&
                                        x.getQuantification_peaks_area() > ms1Trail.getQuantification_peaks_area()
                        ).count();
                    }
                }
            }
            ms1Trail.setdPeakAreaLocalRankInTimeSpan(Math.max(Math.log(20.0 / iPeakAreaRank), 0));

        });
        pbMS1.stop();

        Map<Double,Set<String>> mapRtPepSet = Collections.synchronizedMap(new HashMap<>());
        for(double drt:arrRts)
        {
            Set<String> pepSet = Collections.synchronizedSet(new HashSet<>());
            mapRtPepSet.put(drt,pepSet);
        }
        //为节省内存，一个一个Window process
//        for(int iw=0;iw<iswc.windows.size();iw++)
        int i=0;
        for(MSTwoTrailSet ms2window:iswc.windows)
        {
//test
            i++;
            if(i>2) {
                break;

            }
            List<MatchPeptideWithMultiCharge> listPeptideAdd     = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());
            List<MatchPeptideWithMultiCharge> listPeptideAddTrue = Collections.synchronizedList(new LinkedList<MatchPeptideWithMultiCharge>());
            AtomicInteger icount = new AtomicInteger();


            ProgressBar pb = new ProgressBar("Progress", arrRts.length);
            pb.start();
            Arrays.stream(arrRts).parallel().forEach(
                    y -> {
                        //get all msonetrail in current rt windows
                        long start = System.currentTimeMillis();
                        pb.step();
                        List<MSOneTrail> lMsoneTrail = mapRTMSoneTrail.get(y);

//                        Set<String> uniquePeptideString = Collections.synchronizedSet(new HashSet<>());
//                        Set<String> uniquePeptideString = mapRtPepSet.get(y);


                        if (lMsoneTrail != null) {
                            lMsoneTrail.stream().forEach(
                                    msonetrail -> {
                                        //是否在当前窗口
                                        if (msonetrail.getMz()<=ms2window.mzHigh && msonetrail.getMz()>=ms2window.mzLow &&
                                                msonetrail.getMass() - edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1PPM * pepsMass[0] >= pepsMass[0]) //filter m1 trail which less than the least peptide mass
                                        {
                                            ///////////////////////////////////////
                                            //for each msonetrailFeature get the mz and corresponding window

                                            LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide = new LinkedList<MatchPeptideWithMultiCharge>();
                                            //ms1feature mass +- ppm
                                            int iPosBegin = Arrays.binarySearch(pepsMass, msonetrail.getMass() - edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1PPM * msonetrail.getMass());//pepsMass[0]);//- Utils.H2OMass);

                                            iPosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(pepsMass.length, iPosBegin);//从大于等于thresholdMS1PPM开始
                                            //从小于它的最后一个开始，并考虑所有的相等的情况

                                            for (int iPos = iPosBegin;
                                                 iPos < peps.size() && pepsMass[iPos] < msonetrail.getMass() + edu.uw.waterlooms.peptideMatch.Utils.thresholdMS1PPM * msonetrail.getMass();//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                 iPos++) {
                                                //begin of window with mz of ms1 feature
                                                // if mass matches, then pep match
//                                                List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(msonetrail.getMz(), 1, true);
//                                                for (MSTwoTrailSet window : arrwindow) {
//                                                    MS1MatchPeptide_top6PredictIntensity(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msonetrail, listMaxiumPeptide, peps.get(iPos), window);
                                                MS1MatchPeptide_top6PredictIntensity(icount, peps, bufferedWriter,  MapMS2Intensity, msonetrail, listMaxiumPeptide, peps.get(iPos), ms2window);
//                                                }
                                            }

//                                            mapMSOneFeaturePeptide.put(msonetrail.getId(), listMaxiumPeptide);

                                            if (listMaxiumPeptide.size() > 0) {
                                                if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputALL) {
                                                    listPeptideAddTrue.addAll(listMaxiumPeptide);

                                                } else if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputMax) {
                                                    listPeptideAddTrue.add(listMaxiumPeptide.get(0));
                                                } else if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputSecondMax) {
                                                    if (listMaxiumPeptide.size() > 1) {
                                                        listPeptideAddTrue.add(listMaxiumPeptide.get(1));
                                                    }
                                                }
                                            }
//                                            listPeptideAddTrue.addAll(listMaxiumPeptide);
                                            ims1.addAndGet(listMaxiumPeptide.size());
                                            if (edu.uw.waterlooms.peptideMatch.Utils.bzCheckPeptideWithAllChargeAtSameRT) {
                                                //检查该MS1相同RT且已经匹配的相同PEPTIDE里面所有电荷的二级谱
                                                int iCurrenMS1Z = msonetrail.getZ();
//                                                Set<String> uniquePeptideString = Collections.synchronizedSet( new HashSet<>());
                                                for (MatchPeptideWithMultiCharge matchPep : listMaxiumPeptide) {
                                                    //debug
                                                    if(matchPep.pep.composition.equals("GGDRVPADIR"))
                                                    {
                                                        System.out.println("ms1rt:"+y+"      "+JSONWriter.valueToString(mapRtPepSet.get(y)));
                                                    }
                                                    //是否重复出现的PEPTIDE，若相同只进来一次
//                                                    if (uniquePeptideString.contains(matchPep.pep.composition)) {
                                                    if (mapRtPepSet.get(y).contains(matchPep.pep.composition)) {
                                                        continue;
                                                    }
                                                    for (int iAllC = 1; iAllC <= edu.uw.waterlooms.peptideMatch.Utils.iMS1Charget; iAllC++) {
                                                        //condition is the z is larger than current z
                                                        if (iAllC == iCurrenMS1Z) {
                                                            continue;//means this
                                                        }
                                                        double dPepMZWithC = edu.uw.waterlooms.peptideMatch.Utils.MassToMz(matchPep.pep.mass, iAllC);
                                                        //search the mz in lMsoneTrail,如果存在（可能存在一个或多个）且比当前小就不search 二级谱，否则search，然后在标记该peptide已经被search，下次就不再search

                                                        if (iAllC < iCurrenMS1Z) {
                                                            //if it is matched, it should told the problem there is no matched at successive process(create a map??)
                                                            //get the window
                                                            boolean isE = false;
                                                            int iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC * (1 - edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            for (int iPos = iMS1PosBegin;
                                                                 iPos < lMsoneTrail.size() && lMsoneTrail.get(iPos).getMz() < dPepMZWithC * (1 + edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM);//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                                 iPos++) {

                                                                isE = true;
                                                                break;

                                                            }
                                                            if (isE) {
                                                                continue;
                                                            }

                                                        } else if (iAllC > iCurrenMS1Z) {
                                                            //(the larger z means mz is less than current mz that has been checked before)
                                                            int iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.binarySearch0(lMsoneTrail, 0, lMsoneTrail.size(), dPepMZWithC * (1 - edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM));//pepsMass[0]);//- Utils.H2OMass);

                                                            iMS1PosBegin = edu.uw.waterlooms.peptideMatch.Utils.getiPosBeginOfNearest(lMsoneTrail.size(), iMS1PosBegin);//从大于等于threshold开始
                                                            //从小于它的最后一个开始，并考虑所有的相等的情况
                                                            if (iMS1PosBegin
                                                                    < lMsoneTrail.size() && lMsoneTrail.get(iMS1PosBegin).getMz() < dPepMZWithC * (1 + edu.uw.waterlooms.peptideMatch.Utils.thresholdPPM))//msonetrail.getMass() + Utils.thresholdMS1PPM * pepsMass[iPos];//+ Utils.H2OMass;
                                                            {

                                                                continue;
                                                            }


                                                        }
//                                                        MSTwoTrailSet window = iswc.FindWindowWithMZ(dPepMZWithC, 1);
                                                        List<MSTwoTrailSet> arrwindow = iswc.FindWindowWithMZ(dPepMZWithC, 1, true);
                                                        for (MSTwoTrailSet window : arrwindow) {
                                                            MSOneTrail msOneTrail = new MSOneTrail(-1, dPepMZWithC, msonetrail.getRt(), iAllC, 0, 0, 0);
                                                            LinkedList<MatchPeptideWithMultiCharge> listPeptide = new LinkedList<MatchPeptideWithMultiCharge>();
//                                                            MS1MatchPeptide(icount, peps, bufferedWriter, iswc, MapMS2Intensity, msOneTrail, listPeptide, matchPep.pep, window);
                                                            MS1MatchPeptide(icount, peps, bufferedWriter,  MapMS2Intensity, msOneTrail, listPeptide, matchPep.pep, window);
                                                            if(matchPep.pep.composition.equals("GGDRVPADIR"))
                                                            {
                                                                System.out.println("ms1rt:"+y+"     GGDRVPADIR ");
                                                            }
                                                            listPeptideAdd.addAll(listPeptide);
                                                        }

                                                        //check other charge
                                                    }
                                                    //check whether it has matched before

//                                                    uniquePeptideString.add(matchPep.pep.composition);
                                                    mapRtPepSet.get(y).add(matchPep.pep.composition);
                                                    if(matchPep.pep.composition.equals("GGDRVPADIR"))
                                                    {
                                                        System.out.println("ms1rt:"+y+"      "+JSONWriter.valueToString(mapRtPepSet.get(y)));
                                                    }
                                                }


                                            }
                                        }

                                    }
                            );
                        }
                    }
            );
            pb.stop();


            System.out.println("add listPeptideAddTrue size:" + listPeptideAddTrue.size());

            System.out.println("add CheckPeptideWithAllChargeAtSameRT size:" + listPeptideAdd.size());
            List<MatchPeptideWithMultiCharge> sorted = listPeptideAdd;
            sorted.addAll(listPeptideAddTrue);

            System.out.println("Original Count :" + ims1);

            System.out.println("after combine CheckPeptideWithAllChargeAtSameRT size:" + sorted.size());


            int iCountForFDR = 1;
            int iCountDecoyCurrent = 0;
            Map<Long, Integer> mapIDCount = new HashMap();


            Collections.sort(sorted, Collections.reverseOrder());


            System.out.println("---------OutputALL Begin---------");
            ProgressBar pboutput = new ProgressBar("OutputALL", sorted.size());
            pboutput.start();
            for (MatchPeptideWithMultiCharge mp : sorted) {
                pboutput.step();
                if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputALL) {
                    if (mapIDCount.get(mp.lMsoneid) == null) {
                        mapIDCount.put(mp.lMsoneid, 1);
                    } else {
                        mapIDCount.put(mp.lMsoneid, mapIDCount.get(mp.lMsoneid) + 1);
                    }
                }
                if (mp.strDecoy.equals("decoy")) {
                    iCountDecoyCurrent++;
                }
                bufferedWriter.write((iCountDecoyCurrent + 0.0) / iCountForFDR + "\t");

                if (edu.uw.waterlooms.peptideMatch.Utils.OutputMode == edu.uw.waterlooms.peptideMatch.Utils.OutputModeEnum.OutputALL) {
                    bufferedWriter.write(mp.toString() + "\t");
                    bufferedWriter.write(mapIDCount.get(mp.lMsoneid) + "\n");

                } else {
                    bufferedWriter.write(mp.toString() + "\n");


                }
                if (edu.uw.waterlooms.peptideMatch.Utils.bzOutputConciseAndPeptideFiles) {

                    concisebufferedWriter.write(mp.getConciseString() + "\n");
                    peptidesForpDeep3bufferedWriter.write(mp.pep.composition + "\t\t" + mp.msOneTrail.getZ() + "\n");
                    peptidesForpPrositbufferedWriter.write(mp.pep.composition + "\t27.5\t" + mp.msOneTrail.getZ() + "\tHCD\n");

                }

                iCountForFDR++;

            }
            //for gc
            sorted=null;
            mapIDCount=null;

            pboutput.stop();
            System.gc();

        }

        System.out.println("---------OutputALL Finish---------");

        if (Utils.bzOutputConciseAndPeptideFiles) {
            concisebufferedWriter.flush();
            peptidesForpDeep3bufferedWriter.flush();
            peptidesForpPrositbufferedWriter.flush();
            concisebufferedWriter.close();
            peptidesForpDeep3bufferedWriter.close();
            peptidesForpPrositbufferedWriter.close();
        }


        bufferedWriter.flush();
        bufferedWriter.close();


        System.out.println("Matched with Candidate Peptides ...");

        fullRunStopWatch.stop();
        System.out.println("Elapsed Time in Minutes: " + fullRunStopWatch.getTime(TimeUnit.MINUTES) + " ...");
        System.exit(0);
    }

    private static void MS1MatchPeptide_top6PredictIntensity(AtomicInteger icount, List<edu.uw.waterlooms.peptideMatch.Peptide> peps, BufferedWriter bufferedWriter,
                                                             Map<String, MS2Intensity> MapMS2Intensity, MSOneTrail msonetrail,
                                                             LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide, edu.uw.waterlooms.peptideMatch.Peptide pepsiPos, MSTwoTrailSet window) {
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
                //TEST_DEBU

                ////为获得最大值，或者是加和值，每一次重复计算
                double score = 0.0;


                MatchPeptideWithMultiCharge matchPep = edu.uw.waterlooms.peptideMatch.DbMatch.PepMatchIonCountCrossRTMultiChargeFragmentLoss_top6PredictIntensity(pepsiPos, mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                        Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail,arrbLossCount,arryLossCount
                        ,arrIntensity,ms2Int);
                if (matchPep != null) {
                    matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错

                    if (ms2Int!=null) {
                        matchPep.dPrositRT = ms2Int.dPredictRT;
                    }
                    matchPep.calculateCombineScore();
                    matchPep.setDms2rt(arrMSTwoRts[iRTPos]);


//                    matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);
                    matchPep.iWindowSize = window.iWindowSize;


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


                    }
                    dPeptideMatchScoreSum += matchPep.dMatch;//将在时间内该Peptide所有的可以匹配的值全部加起来
                    if(Utils.PeptideMatchMode==Utils.PeptideMatchEnum.COMBINEMODE) {
                        if (matchPeptdideAll == null) {
                            matchPeptdideAll = matchPep;

                        } else {
                            matchPeptdideAll.combineResult(matchPep);
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
                }
            }
            MatchPeptideWithMultiCharge mp = null;
            double dPm = 0.0;
            if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
            {
                mp = matchPeptdideMax;
                dPm = dPeptideMatchMax;
                if(mp!=null)
                {
                    if(Utils.bzCheckMS2TrailContainMS1Trail) {
                        mp.setStrMS1trailInMS2TrailCount("" + iMS1TrailInMS2TrailCount);
                    }
                }

            } else if(Utils.PeptideMatchMode==Utils.PeptideMatchEnum.COMBINEMODE) {
                mp = matchPeptdideAll;
                dPm = dPeptideMatchScoreSum;
            }


            //peptide search mode: the more ion covered, the more score have
            if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.CoverMax) {

                if (mp != null) {
                    int length = mp.lBIonCovered.stream().distinct().toArray().length +
                            mp.lYIonCovered.stream().distinct().toArray().length;
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
    private static void MS1MatchPeptide(AtomicInteger icount, List<edu.uw.waterlooms.peptideMatch.Peptide> peps, BufferedWriter bufferedWriter,
                                        Map<String, MS2Intensity> MapMS2Intensity, MSOneTrail msonetrail,
                                        LinkedList<MatchPeptideWithMultiCharge> listMaxiumPeptide, edu.uw.waterlooms.peptideMatch.Peptide pepsiPos, MSTwoTrailSet window) {
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




            //get twotrail RT distance between Math.abs(y-v) <= Utils.thresholdRT
            //0725add rt[getStratRt,getEndRt]
            int iTwoRTPosBegin = Arrays.binarySearch(arrMSTwoRts, msonetrail.getStratRt() - Utils.thresholdRT);
            iTwoRTPosBegin = Utils.getiPosBeginOfNearest(arrMSTwoRts.length, iTwoRTPosBegin);
            if (iTwoRTPosBegin > 0) {
                iTwoRTPosBegin = iTwoRTPosBegin - 1;//move to less than position
            }

            for (int iRTPos = iTwoRTPosBegin;
                 iRTPos < arrMSTwoRts.length && arrMSTwoRts[iRTPos] < msonetrail.getEndRt() + Utils.thresholdRT;
                 iRTPos++) {

                ////为获得最大值，或者是加和值，每一次重复计算
                double score = 0.0;

                pepsiPos.GenerateIons();
                double[][] arrIntensity = null;
                MS2Intensity ms2Int = MapMS2Intensity.get(pepsiPos.composition+ msonetrail.getZ());
                if (ms2Int!=null) {
                    arrIntensity = ms2Int.arrdIntensity;
                }



                MatchPeptideWithMultiCharge matchPep = edu.uw.waterlooms.peptideMatch.DbMatch.PepMatchIonCountCrossRTMultiChargeFragmentLoss(pepsiPos, mapRTMStwoTrail.get(arrMSTwoRts[iRTPos]), window.mapRTMSTwoMaxPeakArea.get(arrMSTwoRts[iRTPos]),
                        Config.scoreMode, bufferedWriter, icount.getAndIncrement(), msonetrail,BionMatchTrail,YionMatchTrail,arrbLossCount,arryLossCount
                        ,arrIntensity);
                if (matchPep != null) {
                    matchPep.setMsOneTrail(msonetrail);//必须先设置，否则后面计算combineScore出错

                    if (ms2Int!=null) {
                        matchPep.dPrositRT = ms2Int.dPredictRT;
                    }
                    matchPep.calculateCombineScore();
                    matchPep.setDms2rt(arrMSTwoRts[iRTPos]);


//                    matchPep.iWindowSize = iswc.FindWindowIndex(msonetrail.getMz(), iswc.windows.size() - 1, 0);
                    matchPep.iWindowSize = window.iWindowSize;


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
                    }
                    dPeptideMatchScoreSum += matchPep.dMatch;//将在时间内该Peptide所有的可以匹配的值全部加起来
                    if (matchPeptdideAll == null) {
                        matchPeptdideAll = matchPep;

                    } else {
                        matchPeptdideAll.combineResult(matchPep);
                    }
                    if (matchPep.dMatch > dPeptideMatchMax) {
                        dPeptideMatchMax = matchPep.dMatch;
                        matchPeptdideMax = matchPep;
                    }
                }

            }
            MatchPeptideWithMultiCharge mp = null;
            double dPm;
            if (Utils.PeptideMatchMode == Utils.PeptideMatchEnum.MAXMODE)//strategy two
            {
                mp = matchPeptdideMax;
                dPm = dPeptideMatchMax;
                if(mp!=null)
                {
                    if(Utils.bzCheckMS2TrailContainMS1Trail) {
                        mp.setStrMS1trailInMS2TrailCount("" + iMS1TrailInMS2TrailCount);
                    }
                }

            } else {
                mp = matchPeptdideAll;
                dPm = dPeptideMatchScoreSum;
            }
            //peptide search mode: the more ion covered, the more score have
            if (Utils.PeptideSearchMode == Utils.PeptideSearchIonModeEnum.CoverMax) {

                if (mp != null) {
                    int length = mp.lBIonCovered.stream().distinct().toArray().length +
                            mp.lYIonCovered.stream().distinct().toArray().length;
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


}