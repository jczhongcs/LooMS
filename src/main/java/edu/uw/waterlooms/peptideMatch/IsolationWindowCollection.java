package edu.uw.waterlooms.peptideMatch;


import me.tongfei.progressbar.ProgressBar;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;

public class IsolationWindowCollection implements Serializable {

    // trails are to be sorted by mass
    public ArrayList<MSTwoTrailSet> windows =  new ArrayList();
    public IsolationWindowCollection(ArrayList<MSTwoTrailSet> win) {
        windows = win;
    }
    // constructs isolation windows in file
    // trailInfile is the spectra.mgf outputted by ms2/MS2FeatureDetection
    public IsolationWindowCollection(String trailInfile
    , TreeMap<Double, List<MSOneTrail>> mapRTMSoneFeature)throws IOException {


        String dataFile = trailInfile;
        FileReader freader;
        try {
            freader = new FileReader(dataFile);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
            MSTwoTrailSet cur = new MSTwoTrailSet();
            String line;

            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) { continue; }
                // if see end, add cur to list, next
                if (line.contains("END")) {
                    this.windows.add(cur);
                    continue;
                }

                if (line.contains("START")) {
                    long start = System.currentTimeMillis();

                    cur = new MSTwoTrailSet();
                    line = br.readLine();
                    String[] range = line.split("\\s+");
                    cur.mzLow = Double.parseDouble(range[0]);
                    cur.mzHigh = Double.parseDouble(range[1]);
                    cur.readFromeTrailFile(Config.spectrumFolder +
//                            "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window" +
                            Utils.strinputMS2FragmentInfo +
                            lineType +
                            "_ms2_trails.tsv",mapRTMSoneFeature);
//                    cur.readFromeTrailFile(Config.spectrumFolder +
//                            "output/" +
//                            lineType +
//                            "_window.tsv",mapRTMSoneFeature);
                    lineType ++;
                    long time = System.currentTimeMillis() - start;
                    System.out.println("------------------"+lineType);

                    System.out.println(time);
                    System.out.println(cur.arrMSTwoTrail.size());



                    continue;
                }

            }
            br.close();
            freader.close();
        }
    }



    public IsolationWindowCollection()
    {

    }

    public void IsolationWindowCollection_paralle(String trailInfile
            , TreeMap<Double, List<MSOneTrail>> mapRTMSoneFeature)throws IOException {


        String dataFile = trailInfile;
        FileReader freader;
        try {
            freader = new FileReader(dataFile);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
            MSTwoTrailSet cur = new MSTwoTrailSet();
            String line;

            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) { continue; }
                // if see end, add cur to list, next
                if (line.contains("END")) {
                    this.windows.add(cur);
                    continue;
                }

                if (line.contains("START")) {
//                    long start = System.currentTimeMillis();

                    cur = new MSTwoTrailSet();
                    line = br.readLine();
                    String[] range = line.split("\\s+");
                    cur.mzLow = Double.parseDouble(range[0]);
                    cur.mzHigh = Double.parseDouble(range[1]);
                    cur.iWindowSize = lineType;
//                    cur.readFromeTrailFile(Config.spectrumFolder +
//                            "output/" +
//                            lineType +
//                            "_window.tsv",mapRTMSoneFeature);
                    lineType ++;
//                    long time = System.currentTimeMillis() - start;
//                    System.out.println("------------------"+lineType);
//
//                    System.out.println(time);
//                    System.out.println(cur.arrMSTwoTrail.size());



                    continue;
                }

            }
            br.close();
            freader.close();
        }
        ProgressBar pb = new ProgressBar("Progress", this.windows.size());
        pb.start();
        this.windows.stream().parallel().forEach(x->
        {
            pb.step();
            try {
                x.readFromeTrailFile(Config.spectrumFolder +
                        Utils.strinputMS2FragmentInfo +
                        x.iWindowSize  +
                        "_ms2_trails.tsv",mapRTMSoneFeature);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
        pb.stop();

        this.windows.stream().forEach(x->
        {
            System.out.println(x.iWindowSize+"Window, Ori size:"+x.iOriCount+" Final size:"+x.iFinalCount);
        });
    }


    public IsolationWindowCollection(String trailInfile) throws IOException {

        String dataFile = trailInfile;
        FileReader freader;
        try {
            freader = new FileReader(dataFile);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
            MSTwoTrailSet cur = new MSTwoTrailSet();
            String line;

            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) { continue; }
                // if see end, add cur to list, next
                if (line.contains("END")) {
                    this.windows.add(cur);
                    continue;
                }

                if (line.contains("START")) {
                    long start = System.currentTimeMillis();

                    cur = new MSTwoTrailSet();
                    line = br.readLine();
                    String[] range = line.split("\\s+");
                    cur.mzLow = Double.parseDouble(range[0]);
                    cur.mzHigh = Double.parseDouble(range[1]);
                    cur.readFromeTrailFile(Config.spectrumFolder +
//                            "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01_ms2_trails/Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML_isolation_window" +
                            Utils.strinputMS2FragmentInfo+
                            lineType +
                            "_ms2_trails.tsv");
//                    cur.readFromeTrailFile(Config.spectrumFolder +
//                            "output/" +
//                            lineType +
//                            "_window.tsv");
                    lineType ++;
                    long time = System.currentTimeMillis() - start;
                    System.out.println("------------------"+lineType);

                    System.out.println(time);
                    System.out.println(cur.arrMSTwoTrail.size());

                    continue;
                }

            }
            br.close();
            freader.close();
        }
    }


    // given a precursor peptide mass, find the correct isolation window within collection)
    public MSTwoTrailSet FindWindowWithMZ(double mz, int charge) {
        int len = windows.size();
        int found = FindWindowIndex(mz, len - 1, 0);
        if (found >= 0) {
            return windows.get(found);
        }
        return null;
    }


    public List<MSTwoTrailSet> FindWindowWithMZ(double mz, int charge,boolean bzCorssWindow) {
        List<MSTwoTrailSet> listRes = new ArrayList<>();
        int len = windows.size();
        int found = FindWindowIndex(mz, len - 1, 0);
        if (found >= 0) {
            listRes.add(windows.get(found));
            if(found>0)
            {
                MSTwoTrailSet ms2TrailSetBefore  = windows.get(found-1);
                if(ms2TrailSetBefore.mzHigh>=mz &&  ms2TrailSetBefore.mzLow<=mz)
                {
                    listRes.add(ms2TrailSetBefore);
                }
            }
            if(found<windows.size()-1)
            {
                MSTwoTrailSet ms2TrailSetAfter  = windows.get(found+1);
                if(ms2TrailSetAfter.mzHigh>=mz &&  ms2TrailSetAfter.mzLow<=mz)
                {
                    listRes.add(ms2TrailSetAfter);
                }
            }
        }
        return listRes;
    }

    // given a precursor peptide mass, find the correct isolation window within collection)
    public MSTwoTrailSet FindWindow(double mass, int charge) {
        double mz = Utils.MassToMz(mass, charge);
        int len = windows.size();
        int found = FindWindowIndex(mz, len - 1, 0);
        if (found >= 0) {
            return windows.get(found);
        }
        return null;
    }

    // given mz, return correct isolation window index
    public int FindWindowIndex(double mz, int r, int l) {
        if (r >= l) {
            int mid = l + (r - l) / 2;
            // If the element is present at the
            // middle itself
            if (windows.get(mid).mzHigh >= mz && windows.get(mid).mzLow <= mz)
                return mid;

            // If element is smaller than mid, then
            // it can only be present in left subarray
            if (windows.get(mid).mzLow > mz)
                return FindWindowIndex(mz, mid - 1, l);

            // Else the element can only be present
            // in right subarray
            if (windows.get(mid).mzHigh < mz)
                return FindWindowIndex(mz, r,mid + 1);
        }

        return -1;
    }


    public int FindWindowIndexNoRecur(double mz, int r, int l) {
        while (r >= l) {
            int mid = l + (r - l) / 2;
            // If the element is present at the
            // middle itself
            if (windows.get(mid).mzHigh >= mz && windows.get(mid).mzLow <= mz)
                return mid;

            // If element is smaller than mid, then
            // it can only be present in left subarray
            if (windows.get(mid).mzLow > mz)
            {
                r=mid-1;

            }
//                return FindWindowIndex(mz, mid - 1, l);

            // Else the element can only be present
            // in right subarray
            if (windows.get(mid).mzHigh < mz)
            {
                l=mid+1;
            }
//                return FindWindowIndex(mz, r,mid + 1);
        }

        return -1;
    }


    @Override
    public String toString() {
        return
                "windows=" + windows ;
    }

    public boolean isSmaeMS1Window(double msoneMZ, double msoneMZ1) {
        int len = windows.size();
        int foundFirst = FindWindowIndexNoRecur(msoneMZ, len - 1, 0);
        int foundSecond = FindWindowIndexNoRecur(msoneMZ1, len - 1, 0);
//        int foundFirst = FindWindowIndex(msoneMZ, len - 1, 0);
//        int foundSecond = FindWindowIndex(msoneMZ1, len - 1, 0);
        return (foundFirst>=0) && (foundFirst==foundSecond);

    }
}
