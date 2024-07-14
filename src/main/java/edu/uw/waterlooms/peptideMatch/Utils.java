package edu.uw.waterlooms.peptideMatch;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.round;

public class Utils {

    public static double thresholdPPM = 10.0/1000000;//15 ppm
    public static double thresholdMS1PPM = 10.0/1000000;//30 ppm
    public static double thresholdMS1QualityScore = 0.0;//ms1 feature Quality score threshold
    public static double thresholdRT = 5.0/60;//5.0/60;//5s
    public static double thresholdMS1RT = 5.0/60;//15.0/60;//5/60;//15.0/60;//5s

    public static double thresholdMS1CheckPeakAreaRT = 15.0/60;//15.0/60;//5/60;//15.0/60;//5s

    public static int thresholdPeptideSizeMin = 5;//peptide lenth
    public static double thresholdMassDelton =20;//the difference between ms1 precursor with peptide mass
    public static final double HMass = 1.007825035D;
    public static final double NMass = 14.003074D;
    public static final double N15Mass = 15.000108898D;
    public static final double OMass= 15.99491463D;
    public static final double SMass = 31.9720707D;
    public static final double H2Mass = 2.01565007D;
    public static final double NHMass = 15.010899035D;
    public static final double NH2Mass = 16.01872407D;
    public static final double H2OMass = 18.0105647D;

    public static double C12_MASS = 12;
    public static double C13_MASS = 13.003355;


//    public static final double H2OMass = 0.0;

    public static final double NH3Mass = 17.026549105D;

    public static final double COMass = 27.99491463D;


    public static IonMatchEnum IonMatchMode = IonMatchEnum.MaxAbundancesAreaInRange;
    public static int iCleavageThreshold = 2;
    public static int thresholdPeptideSizeMax = 98;
    public static boolean bzGenerateSemiTrypsin = false;//true;//
    public static boolean bzGenerateMissCleavage = true;//false;

    public static boolean bRemovePeptideFromBothTargetAndDecoy = true;//移除PEPTIDE来源于target and decoy sequence

    public static boolean bzRemoveTrailFromSamePrecursor = true;//false;//remove trail from ms1 precursor
    public static boolean bzJustOneTrailCanMatchAPeptide = false;//true;//
    public static boolean bzCheckMS2TrailContainMS1Trail = true;//false;
    public static int iMS2Charge = 2;//ms2 charge count;

    public static int iMS1Charget = 8;//ms1 charge count;
    public static boolean bzOutputConciseAndPeptideFiles = false;//true;//false;//true;

    public static boolean bzMatchOneRT = false;//true;//false;//true;

    public static boolean bzCalculateThePearsonCorresionWithPredictMSMS = true;//false;//true;//false;//true;
    public static String strpDeepOutfile = "pDeep3/newVersion/pDeep3/result0609.txt";//"pDeep3/pdeep3ResultOutput0502.txt";//"pDeep3/pdeep3ResultWithTargetoutput.txt";
    public static double NoneMatchedIntensity = -10.0;
    public static String strAutoRTOutfile =  "dfUniqueAutoRT.txt";
    public static int iChargeRange = 1;//Create the charge count
    public static int thresholdFragmentCount = 3;//至少要4个以上的fragment匹配
    public static double thresholdPredictRT = 6.5;
    public static boolean bzInputConciseMSonePrecursor = true;
    public static boolean bzCheckMS2FragmentIsotope = true;//检查左侧峰，判定该峰是否是isotope峰
    public static boolean bzCheckMS2FragmentIsotopeCharge = true;//检查右侧峰，判定该峰是否是正确的isotope电荷
    public static boolean bzCheckMS2FragmentLossType = true;//检查丢失峰，如H2O, NH3
    public static int iMS2FragmentTypeCount = 1;//3;//检查fragment峰数量，为3时是检查丢失峰，如H2O, NH3，为1时不检查loss峰
    public static int thresholdMS1ScanSpan = 1;//ms1 combine时，若为2时可以跨一个scan，即2630和2632可以组成一个新的precursor，为1时不能跨scan，为0时则不做combine
    public static double MS1RTSpan = 0.0411;//ms1 rt span
    public static boolean bzOnlyRawdata = false;//true;//只有rawdata  \\\ false;//带xiangyuan data
    public static double thresholdMS1IntensityPercentagofMax =  0.5;//用于precursor合并时，判断下降后上升的丰度值是否大于以前的高峰值的一倍，若是则生成新的precursor
    public static boolean bzCheckPeptideWithAllChargeAtSameRT = true;//用于是否检测相同保留时间是否存在相同的Peptide
    public static boolean bTrypsinRule = false;
    public static boolean bTrypsinRuleWithoutP= true;

    //////////
    //human
    public static String strIsolationWindowCollectionFileName = "isolationWindowRanges1.out";

//    public static String[] strPrositResultFileName=
//            {
//                    "prositHumanR01Result10281t3d.txt"
////                   "prositResult0609.txt"
//            };

    public static String strinputPrecursor ="Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_precursorsNew.tsv";//"Fig4_mouse_cerebellum_MHRM_R01_T0_precursors.tsv";//
    public static String strinputMS2FragmentInfo ="Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window";//"Fig4_mouse_cerebellum_MHRM_R01_T0mzXML/Fig4_mouse_cerebellum_MHRM_R01_T0.mzXML_isolation_window";//
    public static String strRawIsotopeFileName ="ms1istopeRemoveHighChargeCover.csv";//来自raw的isotope"ms1istopeRemoveHighChargeCoverMouseR01.csv";//
    public static String inputFastaFileName="concat-homo_sapines_1target3decoy.fasta";
//    public static String inputFastaFileName="concat-MQ-e1k2DBD-homo_sapiens.fasta";
    public static String strOutputFileName="newResult1028HumanR01.csv";


    /////////////
    //mouse
//    public static String strinputPrecursor ="Fig4_mouse_cerebellum_MHRM_R02_T0_precursors.tsv";//"Fig4_mouse_cerebellum_MHRM_R01_T0_precursors.tsv";//
//    public static String strinputMS2FragmentInfo ="Fig4_mouse_cerebellum_MHRM_R02_T0mzXMLbac/Fig4_mouse_cerebellum_MHRM_R02_T0.mzXML_isolation_window";//"Fig4_mouse_cerebellum_MHRM_R01_T0mzXML/Fig4_mouse_cerebellum_MHRM_R01_T0.mzXML_isolation_window";//
//    public static String strRawIsotopeFileName ="ms1istopeRemoveHighChargeCoverMouseR02.csv";//来自raw的isotope"ms1istopeRemoveHighChargeCoverMouseR01.csv";//
//    public static String inputFastaFileName="concat-MQ-e1k2DBD-uniprot-compressed_true_download_true_format_fasta_query__28_28taxon-2022.08.24-16.43.38.44_1targetVS3Decoy.fasta";
////    public static String inputFastaFileName="concat-MQ-e1k2DBD-uniprot-compressed_true_download_true_format_fasta_query__28_28taxon-2022.08.24-16.43.38.44.fasta";
//    public static String strOutputFileName="newResult1017MouseR2OriTrailNewFeatureRemoveIL1t3d.csv";
//    public static String strIsolationWindowCollectionFileName = "Fig4_mouse_cerebellum_MHRM_R02_T0mzXMLbac/isolationWindowRanges.out";
////    public static String strIsolationWindowCollectionFileName = "isolationWindowRanges.out1";
    public static String[] strPrositResultFileName=
        {
                "prositResultR02Mouse09021_781894.txt",
                "prositResultR02Mouse0902781895_1563788.txt",
                "prositResultR02Mouse09021563789_2345681.txt",
                "prositResultR02Mouse09022345682_3127574.txt",
                "prositResultR02Mouse09023127575_3909468.txt",
                "prositResultR02Mouse09023909469_4691361.txt",
                "prositResultR02Mouse09024691362_5473254.txt",
                "prositResultR02Mouse09025473255_6255147_new.txt",
                "prositResultR02Mouse09026255148_7037041_new.txt",
                "prositResultR02Mouse09027037042_7818933.txt"
                ,"prositMouseR02Result10161t3d.txt"


//                    "prositResult0609.txt"
        };


    //    public static String resultTestPSMFileName= "resultMouseR02_3TrainPSMFileName0922NewModel.csv";//"resultMouseR02TestPSMFileName0903.csv";

    public static boolean bzOutputSharePeaksMax=false;//不输出sharepeakmax

    public static String resultTestPSMFileName=
            "resultsMouseR02test2_PSMFile_24.csv," +
            "resultsMouseR02Train2_PSMFile_24.csv";
//            "resultsHumanR01_combinedfTrain2T0_PSMFile_16.csv,"+
//            "resultsHumanR01_combinedfTrain2T1_PSMFile_16.csv,"+
//            "resultsHumanR01_combinedfTrain2T2_PSMFile_16.csv,"+
//            "resultsHumanR01_combinedfTrain2T3_PSMFile_16.csv,"+
//            "resultsHumanR01_combinedftest2T0_PSMFile_16.csv,"+
//            "resultsHumanR01_combinedftest2T1_PSMFile_16.csv,"+
//            "resultsHumanR01_combinedftest2T2_PSMFile_16.csv,"+
//            "resultsHumanR01_combinedftest2T3_PSMFile_16.csv"
//        "resultsMouseR02_combinedfTrain2T0_PSMFile_07.csv," +
//        "resultsMouseR02_combinedfTrain2T1_PSMFile_07.csv," +
//        "resultsMouseR02_combinedfTrain2T2_PSMFile_07.csv," +
//        "resultsMouseR02_combinedfTrain2T3_PSMFile_07.csv," +
//        "resultsMouseR02_combinedftest2T0_PSMFile.csv_07," +
//        "resultsMouseR02_combinedftest2T1_PSMFile.csv_07," +
//        "resultsMouseR02_combinedftest2T2_PSMFile.csv_07," +
//        "resultsMouseR02_combinedftest2T3_PSMFile.csv_07" ;
        ;//"resultMouseR02TestPSMFileName0903.csv";
//    public static String strOutputtopSharePeaksFileName=  "resultMouseR02_3TrainPSMFileNameNewModel0922.csv";

    public static String strOutputtopSharePeaksFileName=
                    "sharePeaksresultsMouseR02test2_PSMFile_24.csv," +
                    "sharePeaksresultsMouseR02Train2_PSMFile_24.csv";
//                    "SharePeaksresultsHumanR01_combinedfTrain2T0_PSMFile_16.csv,"+
//                    "SharePeaksresultsHumanR01_combinedfTrain2T1_PSMFile_16.csv,"+
//                    "SharePeaksresultsHumanR01_combinedfTrain2T2_PSMFile_16.csv,"+
//                    "SharePeaksresultsHumanR01_combinedfTrain2T3_PSMFile_16.csv,"+
//                    "SharePeaksresultsHumanR01_combinedftest2T0_PSMFile_16.csv,"+
//                    "SharePeaksresultsHumanR01_combinedftest2T1_PSMFile_16.csv,"+
//                    "SharePeaksresultsHumanR01_combinedftest2T2_PSMFile_16.csv,"+
//                    "SharePeaksresultsHumanR01_combinedftest2T3_PSMFile_16.csv";
    //                    "SharePeaksMouseR02_combinedfTrain2T0_PSMFile_07.csv," +
//                    "SharePeaksMouseR02_combinedfTrain2T1_PSMFile_07.csv," +
//                    "SharePeaksMouseR02_combinedfTrain2T2_PSMFile_07.csv," +
//                    "SharePeaksMouseR02_combinedfTrain2T3_PSMFile_07.csv," +
//                    "SharePeaksMouseR02_combinedftest2T0_PSMFile.csv_07," +
//                    "SharePeaksMouseR02_combinedftest2T1_PSMFile.csv_07," +
//                    "SharePeaksMouseR02_combinedftest2T2_PSMFile.csv_07," +
//                    "SharePeaksMouseR02_combinedftest2T3_PSMFile.csv_07" ;
    public static int iSkipScanMS1 =2;//相邻的TRAIL要跨多少个scan

    public static double dThresholdMS2TrailPeakareaForOnePeak = 2.5;//用于保留只有单峰的，且具有比较高丰度的MS2 TRAIL
    public static int iNumPredictIntensityTopN =6;//预测峰前六个
    public static int iDebugoutput =1;//预测峰前六个输出调试

    public static FileWriter fileWriterDubegInfo;
    public static double dThresholdVirturalPrecursorPSMMatchScore = 1.0;//虚拟precursor的psm得分阈值
    public static double dThresholdVirturalPrecursorPSMCoscinSimilarityScore = 0.4;
    public static int iCount1TopPSM =60000;//用于做autort校准

    static {
        try {
            fileWriterDubegInfo = new FileWriter("debugInfo.txt");
//            fileWriterDubegInfo = new FileWriter(Config.fastaFolder+"debugInfo.txt");
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static BufferedWriter bufferedWriterDebugInfo = new BufferedWriter(fileWriterDubegInfo);

    public enum IonMatchEnum {
        FirstInIonset,
        NearestInIon,
        InRangeWithThreshold,
        MaxAbundancesInRange,
        MaxAbundancesAreaInRange

    }

    public enum FragmentLossType {
        H2O,
        NH3,
        none
    }

    public static MutiIonMatchEnum MutiIonMatchMode = MutiIonMatchEnum.None;

    public enum MutiIonMatchEnum {
        Consective,//连续匹配
        BackLargeMassSearch,//假设后面的离子匹配概率低一些//
        ConbineConsectiveAndBackLargeMassSearch,
        None

    }

    public static boolean bzInputSameTargetMS1Precursor = false;// true;// false;//true;
    public static FastaReadEnum FastaReadEnumMode = FastaReadEnum.None;

    public enum FastaReadEnum {
        TrageMode,//连续匹配
        DecoyMode,//假设后面的离子匹配概率低一些//
        None

    }

    public static PeptideMatchEnum PeptideMatchMode = PeptideMatchEnum.MAXMODE;
    //     public static PeptideMatchEnum PeptideMatchMode =   Utils.PeptideMatchEnum.MAXMODE;//strategy one
    public enum PeptideMatchEnum {
        COMBINEMODE,
        MAXMODE
    }

    public static PeptideSearchIonModeEnum PeptideSearchMode = PeptideSearchIonModeEnum.None;

    public enum PeptideSearchIonModeEnum {
        Repeat,
        CoverMax,
        RepeatPlusCoverMax,
        ConsecutiveRepeat,
        None
    }

    public static ScoreSumModeEnum ScoreSumMode = ScoreSumModeEnum.Abundances;

    public enum ScoreSumModeEnum {
        Abundances,
        MassDistance,
        CombineAbundanceAndMassDistance
    }


    public static OutputModeEnum OutputMode = OutputModeEnum.OutputALL;//OutputALL;//OutputMax;

    public enum OutputModeEnum {
        OutputMax,
        OutputSecondMax,
//        OutputSecondMax,
        OutputALL

    }



    public static int iCandidateListSize =10;//30000;// 20;// 30000;



    public static String OutputParametersString() {
        return "Parameter:{" +
                "thresholdPPM=" + thresholdPPM +
                ", thresholdMS1PPM=" + thresholdMS1PPM +
                ", thresholdMS1QualityScore=" + thresholdMS1QualityScore +
                ", thresholdRT=" + thresholdRT +
                ", thresholdMS1RT=" + thresholdMS1RT +
                ", thresholdPeptideSizeMin=" + thresholdPeptideSizeMin +
                ", thresholdMassDelton=" + thresholdMassDelton +
//                ", HMass=" + HMass +
//                ", NMass=" + NMass +
//                ", N15Mass=" + N15Mass +
//                ", OMass=" + OMass +
//                ", SMass=" + SMass +
//                ", H2Mass=" + H2Mass +
//                ", NHMass=" + NHMass +
//                ", NH2Mass=" + NH2Mass +
//                ", H2OMass=" + H2OMass +
//                ", NH3Mass=" + NH3Mass +
//                ", COMass=" + COMass +
                ", IonMatchMode=" + IonMatchMode +
                ", iCleavageThreshold=" + iCleavageThreshold +
                ", thresholdPeptideSizeMax=" + thresholdPeptideSizeMax +
                ", bzGenerateSemiTrypsin=" + bzGenerateSemiTrypsin +
                ", bzGenerateMissCleavage=" + bzGenerateMissCleavage +
                ", bzRemoveTrailFromSamePrecursor=" + bzRemoveTrailFromSamePrecursor +
                ", bzJustOneTrailCanMatchAPeptide=" + bzJustOneTrailCanMatchAPeptide +
                ", bzCheckMS2TrailContainMS1Trail=" + bzCheckMS2TrailContainMS1Trail +
                ", iMS2Charge=" + iMS2Charge +
                ", bzOutputConciseAndPeptideFiles=" + bzOutputConciseAndPeptideFiles +
                ", bzMatchOneRT=" + bzMatchOneRT +
                ", bzCalculateThePearsonCorresionWithPredictMSMS=" + bzCalculateThePearsonCorresionWithPredictMSMS +
                ", strpDeepOutfile='" + strpDeepOutfile + '\'' +
                ", NoneMatchedIntensity=" + NoneMatchedIntensity +
                ", strAutoRTOutfile='" + strAutoRTOutfile + '\'' +
                ", iChargeRange=" + iChargeRange +
                ", MutiIonMatchMode=" + MutiIonMatchMode +
                ", bzInputSameTargetMS1Precursor=" + bzInputSameTargetMS1Precursor +
                ", bzInputConciseMSonePrecursor=" + bzInputConciseMSonePrecursor +
                ", FastaReadEnumMode=" + FastaReadEnumMode +
                ", PeptideMatchMode=" + PeptideMatchMode +
                ", PeptideSearchMode=" + PeptideSearchMode +
                ", ScoreSumMode=" + ScoreSumMode +
                ", OutputMode=" + OutputMode +
                ", iCandidateListSize=" + iCandidateListSize +
                ", thresholdFragmentCount=" + thresholdFragmentCount +
                ", thresholdPredictRT=" + thresholdPredictRT +
                ", bzCheckMS2FragmentIsotope=" + bzCheckMS2FragmentIsotope +
                ", bzCheckMS2FragmentIsotopeCharge=" + bzCheckMS2FragmentIsotopeCharge +
                ", bzCheckMS2FragmentLossType=" + bzCheckMS2FragmentLossType +
                ", iMS2FragmentTypeCount=" + iMS2FragmentTypeCount +
                ", thresholdMS1ScanSpan=" + thresholdMS1ScanSpan +
                ", bzOnlyRawdata=" + bzOnlyRawdata +
                ", thresholdMS1IntensityPercentagofMax=" + thresholdMS1IntensityPercentagofMax +
                ", bzCheckAllChargeWithSameRT=" + bzCheckPeptideWithAllChargeAtSameRT +
                ", iMS1Charget=" + iMS1Charget +
                ", thresholdMS1CheckPeakAreaRT=" + thresholdMS1CheckPeakAreaRT +
                ", bRemovePeptideFromBothTargetAndDecoy=" + bRemovePeptideFromBothTargetAndDecoy +
                ", bTrypsinRule=" + bTrypsinRule +
                ", bTrypsinRuleWithoutP=" + bTrypsinRuleWithoutP +
                ", strinputPrecursor=" + strinputPrecursor +
                ", strinputMS2FragmentInfo=" + strinputMS2FragmentInfo +
                ", strRawIsotopeFileName=" + strRawIsotopeFileName +
                ", inputFastaFileName=" + inputFastaFileName +
                ", strOutputFileName=" + strOutputFileName +
                ", resultTestPSMFileName=" + resultTestPSMFileName +
                ", strOutputtopSharePeaksFileName=" + strOutputtopSharePeaksFileName +
                ", iSkipScanMS1=" + iSkipScanMS1 +
                ", strIsolationWindowCollectionFileName=" + strIsolationWindowCollectionFileName +
                ", strPrositResultFileName=" + strPrositResultFileName +
                ", dThresholdMS2TrailPeakareaForOnePeak=" + dThresholdMS2TrailPeakareaForOnePeak +
                ", bzOutputSharePeaksMax=" + bzOutputSharePeaksMax +
                ", iNumPredictIntensityTopN=" + iNumPredictIntensityTopN +
                ", dThresholdVirturalPrecursorPSMMatchScore=" + dThresholdVirturalPrecursorPSMMatchScore +
                ", dThresholdVirturalPrecursorPSMCoscinSimilarityScore=" + dThresholdVirturalPrecursorPSMCoscinSimilarityScore +
                ", iCount1TopPSM=" + iCount1TopPSM +

                '}';
    }

    public static int getMaxPos(double[] arr) {
        return IntStream.range(0, arr.length).reduce((i, j) -> arr[i] < arr[j] ? j : i).getAsInt();
    }

    public static int getMinPos(double[] arr) {
        return IntStream.range(0, arr.length).reduce((i, j) -> arr[i] > arr[j] ? j : i).getAsInt();
    }


    public static int DeterminePepCharge(double mass) {
        if (mass <= 800) {
            return 2;
        }
        return 1;
    }

    // x is how much you want to round to (1/3 etc)
    public static double RoundToX(double o, double x) {
        if (x == 0) {
            return o;
        }
        double res = round(o / x) * x;
        return res;
    }


    public static double MassToMz(double mass, int z) {
        return mass/z + AminoAcid.ProtonMass;
    }

    public static double MzToMass(double mz, int z) {
        return (mz - AminoAcid.ProtonMass) * z;
    }


    public static int getiPosBeginOfNearest(int length, int iPosBegin) {
        if (iPosBegin < 0)
        {
            int indexOfNearest = ~iPosBegin;

            if (indexOfNearest == length)
            {
                //from time is larger than all elements
                iPosBegin = indexOfNearest;
            }
            else if (indexOfNearest == 0)
            {
                // from time is less than first item
                iPosBegin = 0;
            }
            else
            {
                // from time is between (indexOfNearest - 1) and indexOfNearest
                iPosBegin = indexOfNearest;
            }
        }
        return iPosBegin;
    }


    public static double coscinSimilarity(double[] ms1rt, double[] ms1ins, double[] ms2rt, double[] ms2ins) {
        int ims1 = 0;
        int ims2 = 0;
        double consinSimilarity = 0.0;
        double dotMultiply = 0.0;
        double normA = 0.0;
        double normB = 0.0;

        while (ims1< ms1rt.length || ims2< ms2rt.length)
        {
            for (; ims1== ms1rt.length || (ims2< ms2rt.length) && (ms2rt[ims2]< ms1rt[ims1]); ims2++)
            {
                if(ims2== ms2rt.length ) {
                    break;
                }
                dotMultiply += ms2ins[ims2] *  0.0;
                normA += Math.pow(ms2ins[ims2], 2);
                normB += 0.0;

//                System.out.println("ms1:ms2---0:"+(ims2+1)+"ms1:ms2---0:"+ ms2rt[ims2]);

            }

            for (; ims2== ms2rt.length
                    || (((ims1+1)< ms1rt.length) && ms2rt[ims2] >= ms1rt[ims1 + 1]); ims1++)
            {
                if (ims1== ms1rt.length) {
                    break;
                }
                dotMultiply += 0.0 * ms1ins[ims1];
                normA += 0.0;
                normB += Math.pow(ms1ins[ims1], 2);
//                System.out.println("ms1:ms2---"+(ims1+1)+":0"+"ms1:ms2---"+ ms1rt[ims1]+":0");

            }
            if(ims1< ms1rt.length && ims2< ms2rt.length)
            {
                dotMultiply += ms2ins[ims2] * ms1ins[ims1];
                normA += Math.pow(ms2ins[ims2], 2);
                normB += Math.pow(ms1ins[ims1], 2);

//                System.out.println("ms1:ms2---"+(ims1+1)+":"+(ims2+1)+"ms1:ms2---"+ ms1rt[ims1]+":"+ ms2rt[ims2]);

                ims2++;
                ims1++;

            }

        }
        if((normA>0.0) && (normB>0.0))
        {
            consinSimilarity = dotMultiply / (Math.sqrt(normA) * Math.sqrt(normB));
        }
        return consinSimilarity;
    }
    public static int binarySearch0(List<MSOneTrail> a, int fromIndex, int toIndex,
                                     double key) {
        int low = fromIndex;
        int high = toIndex - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            double midVal = a.get(mid).getMz();

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

    public static String strHtmlTemplatePart1Title = "<!--\n" +
            "    THIS EXAMPLE WAS DOWNLOADED FROM https://echarts.apache.org/examples/zh/editor.html?c=grid-multiple\n" +
            "-->\n" +
            "<!DOCTYPE html>\n" +
            "<html style='height: 100%'>\n" +
            "\n" +
            "<head>\n" +
            "    <meta charset='utf-8'>\n" +
            "</head>\n" +
            "\n" +
            "<body style='height: 100%; margin: 0'>\n" +
            "    <div id='container' style='height: 100%'></div>\n" +
            "    <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/dist/echarts.min.js'></script>\n" +
            "    <!-- Uncomment this line if you want to dataTool extension\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/dist/extension/dataTool.min.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment this line if you want to use gl extension\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts-gl@2/dist/echarts-gl.min.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment this line if you want to echarts-stat extension\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts-stat@latest/dist/ecStat.min.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment this line if you want to use map\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/map/js/china.js'></script>\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/map/js/world.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment these two lines if you want to use bmap extension\n" +
            "        <script type='text/javascript' src='https://api.map.baidu.com/api?v=2.0&ak=<Your Key Here>'></script>\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@{{version}}/dist/extension/bmap.min.js'></script>\n" +
            "        -->\n" +
            "    <script type='text/javascript'>\n" +
            "    var dom = document.getElementById('container');\n" +
            "    var myChart = echarts.init(dom);\n" +
            "    var app = {};\n" +
            "\n" +
            "    var option;\n" +
            "    var differTwoIons = 0.0;\n" +
            "\n" +
            "\n" +
            "\n" +
            "    // prettier-ignore\n" +
            "    // let timeData = [\n" +
            "    //     '2009/10/18 8:00'\n" +
            "    // ];\n" +
            "    // timeData = timeData.map(function (str) {\n" +
            "    //   return str.replace('2009/', '');\n" +
            "    // });\n" +
            "    option = {\n" +
            "        title: {\n" +
            "            text: 'BYIons vs MS2Trail";
    public static String strHtmlTemplatePart1Axi1Max = "',\n" +
            "            left: 'center'\n" +
            "        },\n" +
            "        tooltip: {\n" +
            "            trigger: 'axis',\n" +
            "            axisPointer: {\n" +
            "                animation: true\n" +
            "            }\n" +
            "        },\n" +
            "        legend: {\n" +
            "            data: ['MS2Trail', 'B', 'B+', 'Y', 'Y+'],\n" +
            "            left: 10\n" +
            "        },\n" +
            "        toolbox: {\n" +
            "            feature: {\n" +
            "                dataZoom: {\n" +
            "                    yAxisIndex: 'none'\n" +
            "                },\n" +
            "                restore: {},\n" +
            "                saveAsImage: {}\n" +
            "            }\n" +
            "        },\n" +
            "        axisPointer: {\n" +
            "            link: [{\n" +
            "                xAxisIndex: 'all'\n" +
            "            }]\n" +
            "        },\n" +
            "        dataZoom: [{\n" +
            "                show: true,\n" +
            "                realtime: true,\n" +
            "                start: 30,\n" +
            "                end: 70,\n" +
            "                xAxisIndex: [0, 1]\n" +
            "            },\n" +
            "            {\n" +
            "                type: 'inside',\n" +
            "                realtime: true,\n" +
            "                start: 30,\n" +
            "                end: 70,\n" +
            "                xAxisIndex: [0, 1]\n" +
            "            }\n" +
            "        ],\n" +
            "        grid: [{\n" +
            "                left: 60,\n" +
            "                right: 50,\n" +
            "                top: '6%',\n" +
            "\n" +
            "                height: '40%'\n" +
            "            },\n" +
            "            {\n" +
            "                left: 60,\n" +
            "                right: 50,\n" +
            "                top: '46%',\n" +
            "                height: '40%'\n" +
            "            }\n" +
            "        ],\n" +
            "        xAxis: [{\n" +
            "                position: 'top',\n" +
            "                type: 'value',\n" +
            "                min: 0,\n" +
            "                max: ";
    public static String strHtmlTemplatePart1Axi2Max=",\n" +
            "                //   type: 'category',\n" +
            "                boundaryGap: false,\n" +
            "                axisPointer: {\n" +
            "                    label: {\n" +
            "                      formatter: function (params) {\n" +
            "                         if ((typeof params.seriesData[0]!== 'undefined') && (typeof params.seriesData[0].data !== 'undefined'))\n" +
            "                        {\n" +
            "                            differTwoIons =  params.seriesData[0].data[0] ;\n" +
            "                        }\n" +
            "\n" +
            "                        return (\n" +
            "                          'mass  ' +\n" +
            "                          params.value  +\n" +
            "                          (((typeof params.seriesData[0]!== 'undefined') && (typeof params.seriesData[0].data !== 'undefined')) ? ' NO.：' + params.seriesData[0].data[2] : '')\n" +
            "                        );\n" +
            "                      }\n" +
            "                    }\n" +
            "                },\n" +
            "                axisLine: { onZero: true } //,\n" +
            "                //   data: timeData\n" +
            "            },\n" +
            "            {\n" +
            "                gridIndex: 1,\n" +
            "                //     type: 'category',\n" +
            "                type: 'value',\n" +
            "                min: 0,\n" +
            "                max:";
    public static String strHtmlTemplatePart2Ms2trailData = ",\n" +
            "                boundaryGap: false,\n" +
            "                axisPointer: {\n" +
            "                    label: {\n" +
            "                      formatter: function (params) {\n" +
            "\n" +
            "                        return (\n" +
            "                          'mass  ' +\n" +
            "                          params.value  +'\\n'+\n" +
            "                          'diff:'+(Math.abs(params.value-differTwoIons))+\n" +
            "                          (((Math.abs(params.value-differTwoIons)*1000000)/params.value)<20 ? ' Matched': 'NoMatch')\n" +
            "                        );\n" +
            "                      }\n" +
            "                    }\n" +
            "                },\n" +
            "                axisLine: { onZero: true } //,\n" +
            "\n" +
            "                //      data: timeData\n" +
            "            }\n" +
            "        ],\n" +
            "        yAxis: [{\n" +
            "                name: 'MS2Trail',\n" +
            "                type: 'value'\n" +
            "                // ,\n" +
            "                // max: 2000\n" +
            "            },\n" +
            "            {\n" +
            "                gridIndex: 1,\n" +
            "                name: 'BYIons',\n" +
            "                type: 'value',\n" +
            "                inverse: true\n" +
            "                 ,\n" +
            "                max: 600\n" +
            "            }\n" +
            "        ],\n" +
            "        series: [\n" +
            "\n" +
            "            {\n" +
            "                name: 'MS2Trail',\n" +
            "                type: 'bar',\n" +
            "                symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#0001ff'\n" +
            "                    }\n" +
            "                },\n" +
            "                // label: {\n" +
            "                //   normal: {\n" +
            "                //     show: true,\n" +
            "                //     color: '#323232',\n" +
            "                //     position: 'top'\n" +
            "                //   }\n" +
            "                // },\n" +
            "                // prettier-ignore\n" +
            "                data:";

    public static String strHtmlTemplatePart3BData ="},\n" +
            "            {\n" +
            "                name: 'B',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                        // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'blue',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#49a1ff'\n" +
            "                    }\n" +
            "                },\n" +
            "                // label: {\n" +
            "                //   normal: {\n" +
            "                //     show: true,\n" +
            "                //     color: '#323232',\n" +
            "                //     position: 'top'\n" +
            "                //   }\n" +
            "                // },\n" +
            "                // prettier-ignore\n" +
            "                data:";
    public static String strHtmlTemplatePart3B2ChargeData ="},\n" +
            "            {\n" +
            "                name: 'B+',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#ff00ff'\n" +
            "                    }\n" +
            "                },\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                        // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'red',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                // prettier-ignore\n" +
            "                data: ";
    public static String strHtmlTemplatePart3YData =" },\n" +
            "\n" +
            "            {\n" +
            "                name: 'Y',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#4ff0a1'\n" +
            "                    }\n" +
            "                },\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                        // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'red',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                // prettier-ignore\n" +
            "\n" +
            "                data:";

    public static String strHtmlTemplatePart3Y2ChargeData ="},\n" +
            "            {\n" +
            "                name: 'Y+',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#ff8f00'\n" +
            "                    }\n" +
            "                },\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                       // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'blue',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                // prettier-ignore\n" +
            "                // \n" +
            "                data: ";
    public static String strHtmlTemplatePart4End ="}\n" +
            "        ]\n" +
            "    };\n" +
            "\n" +
            "    if (option && typeof option === 'object') {\n" +
            "        myChart.setOption(option);\n" +
            "    }\nfunction isMatched(dataMass)\n" +
            "    {\n" +
            "//        console.log(option.series[0].data);\n" +
            "        var dinteralDiffThreshold = dataMass * 10/1000000;\n" +
            "        var bMatched  = iterativeFunction(option.series[0].data,dataMass,dinteralDiffThreshold);\n" +
            "\n" +
            "        return bMatched;\n" +
            "\n" +
            "    }\n" +
            "    function iterativeFunction (arr, x, interalvalue) {\n" +
            "  \n" +
            "        var start=0, end=arr.length-1;\n" +
            "\n" +
            "             \n" +
            "        // Iterate while start not meets end\n" +
            "        while (start<=end){\n" +
            "             \n" +
            "\n" +
            "            // Find the mid index\n" +
            "            var mid=Math.floor((start + end)/2);\n" +
            "      \n" +
            "            // serial data is 3 array, If element is present at mid, return True\n" +
            "            // if (arr[0][mid]===x) return mid;\n" +
            "            if (Math.abs(arr[mid][0]-x)<interalvalue) return mid;\n" +
            "            // Else look in left or right half accordingly\n" +
            "            else if (arr[mid][0] < x)\n" +
            "                 start = mid + 1;\n" +
            "            else\n" +
            "                 end = mid - 1;\n" +
            "        }\n" +
            "      \n" +
            "        return -start;\n" +
            "    }\n" +
            "    function formatlab(d)\n" +
            "    {\n" +
            "        var match_h2o = isMatched(d.data[0]-18.0105647);\n" +
            "        var match_h = isMatched(d.data[0]-1.007825035);\n" +
            "        var match_nh3 = isMatched(d.data[0]-17.026549105);\n" +
            "\n" +
            "        return '{a|'+d.seriesName +' \\n'+d.data[2]+'}{x|'\n" +
            "        +((isMatched(d.data[0])>0)?'\\nMatch':'')\n" +
            "        +((match_h2o>0)?('\\nM-H2O pos:'+(match_h2o+1)):'')\n" +
            "        +((match_h>0)?('\\nM-H pos:'+(match_h+1)):'')\n" +
            "        +((match_nh3>0)?('\\nM-NH3 pos:'+(match_nh3+1)):'')\n" +
            "        +'}';\n" +
            "    }" +

            "    </script>\n" +
            "</body>\n" +
            "\n" +
            "</html>";

}
