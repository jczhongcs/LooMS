package edu.uw.waterlooms;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.ms1.FeatureDetect;
import edu.uw.waterlooms.msutil.OpenMzxml;
import edu.uw.waterlooms.peptideMatch.IsolationWindowCollection;
import edu.uw.waterlooms.peptideMatch.MSOneTrail;
import edu.uw.waterlooms.peptideMatch.MSTwoTrailSet;
import me.tongfei.progressbar.ProgressBar;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class GenerateMS2Trail {
    public Map<Integer, ArrayList<XIC>> mapWindowXICs;
    IsolationWindowCollection iswc = new IsolationWindowCollection();

    public GenerateMS2Trail()
    {
        mapWindowXICs = Collections.synchronizedMap(new HashMap<>());

    }


    public void featureDetectShouldReturnAnArrayListOfXICs(String dataDirectory,OpenMzxml openMzxml,
                                                           TreeMap<Double, List<MSOneTrail>> mapRTMSoneFeature) throws IOException {
        long startTime = System.currentTimeMillis();
        /* Given */

//        String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
//        String dataDirectory = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/" ;
//        String mzXMLFileName = "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML";
//        String mzxmlFile = dataDirectory + mzXMLFileName;
        int rtNextPeakTolSec = 10;  // 5sec
        double mzTolerancePPMStrict = 10e-6; // 10ppm
        int MIN_PEAKNUM_Strict = 2;
        double rtMaxRangeSec = 30; // 30sec
        double intensityNextPeakPercentageTol = 0.9; // 90%
        boolean considerLessConfidentTrails = true;
//        String pathStr = dataDirectory+mzXMLFileName.replaceAll("\\.","")+"/";
////        System.out.println(pathStr);
//        final Path directory = Paths.get(pathStr);
//
//        if (Files.notExists(directory)){
//            Files.createDirectories(directory);
//            System.out.println("Success Create the directory"+pathStr);
//        }
        /* When */
//        OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
        FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS2);

//        FileWriter fileWriterTest = new FileWriter(pathStr + "isolationWindowRanges.out");
        FileWriter fileWriterTest = new FileWriter(dataDirectory + "isolationWindowRanges.out");
        BufferedWriter bufferedWriterTest = new BufferedWriter(fileWriterTest);

        /* Then */
        featureDetect.init_MS2();
//        ArrayList<IsolationWindow> isolationWindowlList = new ArrayList<>();

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

//                    isolationWindowlList.add(new IsolationWindow(dMZlo,dMZhi,isolationWindowlList.size()));
                    MSTwoTrailSet cur = new MSTwoTrailSet();
                    cur.mzLow = dMZlo;
                    cur.mzHigh = dMZhi;
                    cur.iWindowSize = iswc.windows.size();
                    iswc.windows.add(cur);

                    bufferedWriterTest.write("START\n"
                            +dMZlo
                            +" "
                            +dMZhi
                            +"\n"
                            +"END\n");
                }
                if (bAllWindows)
                    break;
                iSkipScanMS1++;
            }
//            iWindowSize = isolationWindowlList.size();
        }

        bufferedWriterTest.flush();
        bufferedWriterTest.close();
//        fileWriterTest.flush();
        fileWriterTest.close();

        if(iswc.windows.size()==0)
            System.out.println("This mzXML file may be not DIA data");

        ProgressBar pb = new ProgressBar("Progress", iswc.windows.size());
        pb.start();
        iswc.windows.stream().parallel().forEach(x->
        {
            pb.step();
            MSTwoTrailSet isolationWindow = x;
            featureDetect.detectMS2Features_method2ForWindow(x, isolationWindow.mzLow, isolationWindow.mzHigh, rtNextPeakTolSec, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeSec, intensityNextPeakPercentageTol, considerLessConfidentTrails,mapRTMSoneFeature);
//            ArrayList<XIC> trailsStrictTol = featureDetect.getTrailsStrictTol();

            //set value to ms2trail


//            ArrayList<XIC> trailsLooseTol = featureDetect.getTrailsLooseTol();
//            try {
////                featureDetect.writeMS2TrailsData(pathStr+mzXMLFileName + "_isolation_window" + x.iNo + "_ms2_trails.tsv", trailsStrictTol,trailsLooseTol);
//
//                featureDetect.writeMS2TrailsData(pathStr+mzXMLFileName + "_isolation_window" + x.iNo + "_ms2_trails.tsv", trailsStrictTol);
//            } catch (IOException e) {
//                throw new RuntimeException(e);
//            }
            System.gc();
        });
        pb.stop();
        long endTime = System.currentTimeMillis();
        System.out.println ("Runtime: " + (endTime - startTime) + "ms");
    }
    public void featureDetectShouldReturnAnArrayListOfXICs(String dataDirectory,String mzXMLFileName) throws IOException {
        long startTime = System.currentTimeMillis();
        /* Given */

//        String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
//        String dataDirectory = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/" ;
//        String mzXMLFileName = "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML";
        String mzxmlFile = dataDirectory + mzXMLFileName;
        int rtNextPeakTolSec = 10;  // 5sec
        double mzTolerancePPMStrict = 10e-6; // 10ppm
        int MIN_PEAKNUM_Strict = 2;
        double rtMaxRangeSec = 30; // 30sec
        double intensityNextPeakPercentageTol = 0.9; // 90%
        boolean considerLessConfidentTrails = true;
        String pathStr = dataDirectory+mzXMLFileName.replaceAll("\\.","")+"/";
//        System.out.println(pathStr);
        final Path directory = Paths.get(pathStr);

        if (Files.notExists(directory)){
            Files.createDirectories(directory);
            System.out.println("Success Create the directory"+pathStr);
        }
        /* When */
        OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
        FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS2);

        FileWriter fileWriterTest = new FileWriter(pathStr + "isolationWindowRanges.out");
        BufferedWriter bufferedWriterTest = new BufferedWriter(fileWriterTest);

        /* Then */
        featureDetect.init_MS2();
        ArrayList<IsolationWindow> isolationWindowlList = new ArrayList<>();
//        if(openMzxml.num2scan!=null)
//        {
//            for (int i:openMzxml.num2scan.keySet()) {
//                for (int scan2num : openMzxml.num2scan.get(i).getChildScans()) {
//
//                    isolationWindowlList.add(new IsolationWindow(openMzxml.num2scan2.get(scan2num).getPrecursor().getMzRange().getLo()
//                            ,openMzxml.num2scan2.get(scan2num).getPrecursor().getMzRange().getHi()));
//
//                    bufferedWriterTest.write("START\n"
//                    +openMzxml.num2scan2.get(scan2num).getPrecursor().getMzRange().getLo()
//                    +" "
//                    +openMzxml.num2scan2.get(scan2num).getPrecursor().getMzRange().getHi()
//                    +"\n"
//                    +"END\n");
//
//
//                }
//                break;
//            }
//        }

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

                    isolationWindowlList.add(new IsolationWindow(dMZlo,dMZhi,isolationWindowlList.size()));


                    bufferedWriterTest.write("START\n"
                            +dMZlo
                            +" "
                            +dMZhi
                            +"\n"
                            +"END\n");
                }
                if (bAllWindows)
                    break;
                iSkipScanMS1++;
            }
            iWindowSize = isolationWindowlList.size();
        }

        bufferedWriterTest.flush();
        bufferedWriterTest.close();
//        fileWriterTest.flush();
        fileWriterTest.close();

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



//        for (int i = 0; i < 2; i++) {
////        for (int i = 0; i < isolationWindowlList.size(); i++) {
//
//            IsolationWindow isolationWindow = isolationWindowlList.get(i);
//            featureDetect.detectMS2Features_method2(mzxmlFile, isolationWindow.mzLow, isolationWindow.mzHigh, rtNextPeakTolSec, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeSec, intensityNextPeakPercentageTol, considerLessConfidentTrails);
//            ArrayList<XIC> trailsStrictTol = featureDetect.getTrailsStrictTol();
//            ArrayList<XIC> trailsLooseTol = featureDetect.getTrailsLooseTol();
//            featureDetect.writeMS2TrailsData(pathStr+mzXMLFileName + "_isolation_window" + i + "_ms2_trails.tsv", trailsStrictTol,trailsLooseTol);
//            System.gc();
//        }
        ProgressBar pb = new ProgressBar("Progress", isolationWindowlList.size());
        pb.start();
        isolationWindowlList.stream().parallel().forEach(x->
        {
            pb.step();
            IsolationWindow isolationWindow = x;
            featureDetect.detectMS2Features_method2(mzxmlFile, isolationWindow.mzLow, isolationWindow.mzHigh, rtNextPeakTolSec, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeSec, intensityNextPeakPercentageTol, considerLessConfidentTrails);
            ArrayList<XIC> trailsStrictTol = featureDetect.getTrailsStrictTol();
//            ArrayList<XIC> trailsLooseTol = featureDetect.getTrailsLooseTol();
            try {
//                featureDetect.writeMS2TrailsData(pathStr+mzXMLFileName + "_isolation_window" + x.iNo + "_ms2_trails.tsv", trailsStrictTol,trailsLooseTol);
                featureDetect.writeMS2TrailsData(pathStr+mzXMLFileName + "_isolation_window" + x.iNo + "_ms2_trails.tsv", trailsStrictTol);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
            System.gc();
        });
        pb.stop();
        long endTime = System.currentTimeMillis();
        System.out.println ("Runtime: " + (endTime - startTime) + "ms");
    }

    //java -Xmx96g -cp  /data/waterlooms-1.0-SNAPSHOT-jar-with-dependencies.jar edu.uw.waterlooms.GenerateMS2Trail /data/mouse/ Fig4_mouse_cerebellum_MHRM_R02_T0.mzXML
    public static void main(String[] args) throws IOException {

        GenerateMS2Trail gen = new GenerateMS2Trail();

        gen.featureDetectShouldReturnAnArrayListOfXICs(args[0],args[1]);//directory and mzXML filename

    }

    public ArrayList<IsolationWindow>   ReadfeatureDetectShouldReturnAnArrayListOfXICs(String dataDirectory,String mzXMLFileName) throws IOException {
        long startTime = System.currentTimeMillis();

        String mzxmlFile = dataDirectory + mzXMLFileName;
        int rtNextPeakTolSec = 10;  // 5sec
        double mzTolerancePPMStrict = 10e-6; // 10ppm
        int MIN_PEAKNUM_Strict = 2;
        double rtMaxRangeSec = 30; // 30sec
        double intensityNextPeakPercentageTol = 0.9; // 90%
        boolean considerLessConfidentTrails = true;
        String pathStr = dataDirectory+mzXMLFileName.replaceAll("\\.","")+"/";
//        System.out.println(pathStr);
        final Path directory = Paths.get(pathStr);

        if (Files.notExists(directory)){
            Files.createDirectories(directory);
            System.out.println("Success Create the directory"+pathStr);
        }
        /* When */
        OpenMzxml openMzxml = new OpenMzxml(mzxmlFile);
        FeatureDetect featureDetect = new FeatureDetect(openMzxml, FeatureDetect.DetectionType.MS2);

        FileWriter fileWriterTest = new FileWriter(pathStr + "isolationWindowRanges.out");
        BufferedWriter bufferedWriterTest = new BufferedWriter(fileWriterTest);

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

                    isolationWindowlList.add(new IsolationWindow(dMZlo,dMZhi,isolationWindowlList.size()));


                    bufferedWriterTest.write("START\n"
                            +dMZlo
                            +" "
                            +dMZhi
                            +"\n"
                            +"END\n");
                }
                if (bAllWindows)
                    break;
                iSkipScanMS1++;
            }
            iWindowSize = isolationWindowlList.size();
        }

        bufferedWriterTest.flush();
        bufferedWriterTest.close();
//        fileWriterTest.flush();
        fileWriterTest.close();

        if(isolationWindowlList.size()==0)
            System.out.println("This mzXML file may be not DIA data");


        ProgressBar pb = new ProgressBar("Progress", isolationWindowlList.size());
        pb.start();
        isolationWindowlList.stream().parallel().forEach(x->
        {
            pb.step();
            IsolationWindow isolationWindow = x;
            featureDetect.detectMS2Features_method2(mzxmlFile, isolationWindow.mzLow, isolationWindow.mzHigh, rtNextPeakTolSec, mzTolerancePPMStrict, MIN_PEAKNUM_Strict, rtMaxRangeSec, intensityNextPeakPercentageTol, considerLessConfidentTrails);
            ArrayList<XIC> trailsStrictTol = featureDetect.getTrailsStrictTol();
//            ArrayList<XIC> trailsLooseTol = featureDetect.getTrailsLooseTol();

            mapWindowXICs.put(x.iNo,trailsStrictTol);
            System.gc();
        });
        pb.stop();
        long endTime = System.currentTimeMillis();
        System.out.println ("Runtime: " + (endTime - startTime) + "ms");
        return  isolationWindowlList;
    }

}
