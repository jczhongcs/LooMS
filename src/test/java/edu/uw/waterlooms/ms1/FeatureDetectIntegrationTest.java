package edu.uw.waterlooms.ms1;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.msutil.OpenMzxml;
import org.junit.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Properties;

public class FeatureDetectIntegrationTest {

//    @Test
    //this CODE generate the ms2 trail info zjc
    public void featureDetectShouldReturnAnArrayListOfXICs() throws IOException {
        long startTime = System.currentTimeMillis();
        /* Given */

//        String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
        String dataDirectory = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/" ;
        String mzXMLFileName = "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML";
        String mzxmlFile = dataDirectory + mzXMLFileName;
        int rtNextPeakTolSec = 10;  // 5sec
        double mzTolerancePPMStrict = 10e-6; // 10ppm
        int MIN_PEAKNUM_Strict = 2;
        double rtMaxRangeSec = 30; // 30sec
        double intensityNextPeakPercentageTol = 0.9; // 90%
        boolean considerLessConfidentTrails = true;

        /* When */
        OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
        FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS2);

        /* Then */
        featureDetect.init_MS2();
        ArrayList<IsolationWindow> isolationWindowlList = new ArrayList<>();
        if(openMzxml.num2scan!=null)
        {
            boolean bfirst = true;
            boolean bAllWindows = false;
            double dStratMZ=0.0;
            int iSkipScanMS1 = 0;
            int iWindowSize = 0;

            for (int i:openMzxml.num2scan.keySet()) {

                for (int scan2num : openMzxml.num2scan.get(i).getChildScans()) {
                    double dMZlo = openMzxml.num2scan2.get(scan2num).getPrecursor().getMzRange().getLo();
                    double dMZhi = openMzxml.num2scan2.get(scan2num).getPrecursor().getMzRange().getHi();
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

                    isolationWindowlList.add(new IsolationWindow(dMZlo,dMZhi));
                }
                if (bAllWindows)
                    break;
                iSkipScanMS1++;
            }
            iWindowSize = isolationWindowlList.size();
        }
        if(isolationWindowlList.size()==0)
            System.out.println("This mzXML file may be not DIA data");

//        isolationWindowlList.add(new IsolationWindow(350,387));
//        isolationWindowlList.add(new IsolationWindow(386,416));
//        isolationWindowlList.add(new IsolationWindow(415,439));
//        isolationWindowlList.add(new IsolationWindow(438,462));
//        isolationWindowlList.add(new IsolationWindow(461,483));
//        isolationWindowlList.add(new IsolationWindow(482,505));
//        isolationWindowlList.add(new IsolationWindow(504,525));
//        isolationWindowlList.add(new IsolationWindow(524,548));
//        isolationWindowlList.add(new IsolationWindow(547,568));
//        isolationWindowlList.add(new IsolationWindow(567,591));
//        isolationWindowlList.add(new IsolationWindow(590,614));
//        isolationWindowlList.add(new IsolationWindow(613,638));
//        isolationWindowlList.add(new IsolationWindow(637,664));
//        isolationWindowlList.add(new IsolationWindow(663,690));
//        isolationWindowlList.add(new IsolationWindow(689,719));
//        isolationWindowlList.add(new IsolationWindow(718,753));
//        isolationWindowlList.add(new IsolationWindow(752,790));
//        isolationWindowlList.add(new IsolationWindow(789,832));
//        isolationWindowlList.add(new IsolationWindow(831,884));
//        isolationWindowlList.add(new IsolationWindow(883,955));
//        isolationWindowlList.add(new IsolationWindow(954,1057));
//        isolationWindowlList.add(new IsolationWindow(1056,1650));

//        isolationWindowlList.add(new IsolationWindow(350,379));
//        isolationWindowlList.add(new IsolationWindow(378,402));
//        isolationWindowlList.add(new IsolationWindow(401,421));
//        isolationWindowlList.add(new IsolationWindow(420,438));
//        isolationWindowlList.add(new IsolationWindow(437,455));
//        isolationWindowlList.add(new IsolationWindow(454,470));
//        isolationWindowlList.add(new IsolationWindow(469,486));
//        isolationWindowlList.add(new IsolationWindow(485,501));
//        isolationWindowlList.add(new IsolationWindow(500,517));
//        isolationWindowlList.add(new IsolationWindow(516,533));
//        isolationWindowlList.add(new IsolationWindow(532,549));
//        isolationWindowlList.add(new IsolationWindow(548,564));
//        isolationWindowlList.add(new IsolationWindow(563,580));
//        isolationWindowlList.add(new IsolationWindow(579,597));
//        isolationWindowlList.add(new IsolationWindow(596,614));
//        isolationWindowlList.add(new IsolationWindow(613,631));
//        isolationWindowlList.add(new IsolationWindow(630,650));
//        isolationWindowlList.add(new IsolationWindow(649,669));
//        isolationWindowlList.add(new IsolationWindow(668,688));
//        isolationWindowlList.add(new IsolationWindow(687,709));
//        isolationWindowlList.add(new IsolationWindow(708,732));
//        isolationWindowlList.add(new IsolationWindow(731,758));
//        isolationWindowlList.add(new IsolationWindow(757,785));
//        isolationWindowlList.add(new IsolationWindow(784,815));
//        isolationWindowlList.add(new IsolationWindow(814,848));
//        isolationWindowlList.add(new IsolationWindow(847,889));
//        isolationWindowlList.add(new IsolationWindow(888,939));
//        isolationWindowlList.add(new IsolationWindow(938,1002));
//        isolationWindowlList.add(new IsolationWindow(1001,1092));
//        isolationWindowlList.add(new IsolationWindow(1091,1650));
        File directory = new File(dataDirectory+mzXMLFileName+"/");
        if (! directory.exists()){
            directory.mkdir();
            // If you require it to make the entire directory path including parents,
            // use directory.mkdirs(); here instead.
        }
        for (int i = 0; i < isolationWindowlList.size(); i++) {
            IsolationWindow isolationWindow = isolationWindowlList.get(i);
            featureDetect.detectMS2Features_method2(mzxmlFile, isolationWindow.mzLow, isolationWindow.mzHigh, rtNextPeakTolSec, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeSec, intensityNextPeakPercentageTol, considerLessConfidentTrails);
            ArrayList<XIC> trailsStrictTol = featureDetect.getTrailsStrictTol();
//            ArrayList<XIC> trailsLooseTol = featureDetect.getTrailsLooseTol();
            featureDetect.writeMS2TrailsData(dataDirectory+mzXMLFileName+"/"+mzXMLFileName + "_isolation_window" + i + "_ms2_trails.tsv", trailsStrictTol);
            System.gc();
        }
        long endTime = System.currentTimeMillis();
        System.out.println ("Runtime: " + (endTime - startTime) + "ms");
    }

//   @Test
public void matchMS1TrailsWithFeatureDetectionResult() throws IOException {
    /* Given */
    String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
       dataDirectory = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/" ;

       String featureDetectParamFile = dataDirectory + "/featuredetect.params";
    String rawFileName = "mouse/Fig4_mouse_cerebellum_MHRM_R02_T0.mzXML";
    String mzxmlFile = dataDirectory + rawFileName;

    /* When */
    OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
    FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS1);

    /* Then */
    featureDetect.detectFeatures(mzxmlFile);
}
//@Test
    public void searchIsotopeWithSignalTest() throws IOException {
        /* Given */
        String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
        String featureDetectParamFile = dataDirectory + "/featuredetect.params";

//            String rawFileName = "toy.mzXML";
            String rawFileName =  "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML";
//            String rawFileName =  "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R03.mzXML";
//            String rawFileName =  "Fig4_mouse_cerebellum_MHRM_R02_T0.mzXML";
            dataDirectory = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/" ;
//            dataDirectory = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/mouse/" ;
        String mzxmlFile = dataDirectory + rawFileName;

        /* When */
        OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
        FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS1);
        featureDetect.init_MS1();
        /* Then */
//        featureDetect.detectFeatures(mzxmlFile);
            featureDetect.searchIsotopeWithSignal(4,8);
    }

    /**
     * Given an Absolute Path to a FeatureDetect Params File, return the Parameters in a Property
     * @param featureDetectParamFile String with absolute value of path
     * @return Properties of parameters file
     */
    private Properties parseFeatureDetectParamsFile(String featureDetectParamFile){
        Properties featureDetectParams = new Properties();
        try {
            FileReader fin = new FileReader(featureDetectParamFile);
            featureDetectParams.load(fin);
        } catch (IOException $e) {
            System.err.printf($e.getMessage());
            System.exit(1);
        }
        return featureDetectParams;
    }

}