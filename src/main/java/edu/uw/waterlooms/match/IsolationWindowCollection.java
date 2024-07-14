/*
package edu.uw.waterlooms.match;
import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;
import java.util.*;
import java.io.*;

public class IsolationWindowCollection {
    // trails are to be sorted by mass
    public ArrayList<IsolationWindow> windows =  new ArrayList();

    // constructs isolation windows in file
    // trailInfile is the spectra.mgf outputted by ms2/MS2FeatureDetection
    public IsolationWindowCollection(String trailInfile) throws IOException{

        String dataFile = trailInfile;
        // create new isolation window
        IsolationWindow cur = new IsolationWindow(0,0);
        XIC curTrail = new XIC(0);
        int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail

        int startCounter = 0;

        FileInputStream inputStream = null;
        Scanner sc = null;
        try {
            inputStream = new FileInputStream(dataFile);
            sc = new Scanner(inputStream, "UTF-8");
            while (sc.hasNextLine()) {
                String line = sc.nextLine();
                if (line.isEmpty()) { continue; }
                if (line.contains("END")){
                    // Sort Trails Here
                    Collections.sort(cur.xics);
                    this.windows.add(cur);
                    continue;
                }
                if (line.contains("START")) {
                    // For Progress
                    System.out.println(startCounter + "/" + "22");
                    startCounter += 1;

                    cur = new IsolationWindow(0, 0);
                    line = sc.nextLine();
                    String[] range = line.split("\\s+");
                    cur.mzLow = Double.parseDouble(range[0]);
                    cur.mzHigh = Double.parseDouble(range[1]);
                    cur.xics = new ArrayList<>();
                    lineType = 0;
                    continue;
                }
                if (lineType == 0) {
                    String[] mzLine = line.split("\\s+");
                    curTrail.setMZAtMaxIntensity(Double.parseDouble(mzLine[0]));
                    curTrail.setRtAtMaxIntensity(Double.parseDouble(mzLine[1]));
                    lineType = 1;
                    continue;
                }
                if (lineType == 1) {
                    String[] linelst = line.split("\\s+");
                    curTrail.setIntensities(
                        Arrays.stream(linelst)
                        .mapToDouble(Double::parseDouble).toList()
                    );
                    lineType = 2;
                    continue;
                }
                if (lineType == 2) {
                    String[] linelst = line.split("\\s+");
                    curTrail.ParseRts(linelst);
                    // add to list
                    cur.PutRtMax(curTrail.rts, curTrail.intensities);
                    cur.xics.add(curTrail);
                    curTrail = new IsotopeTrail();
                    lineType = 0;
                    continue;
                }
                // System.out.println(line);
            }
            // note that Scanner suppresses exceptions
            if (sc.ioException() != null) {
                throw sc.ioException();
            }
        } finally {
            if (inputStream != null) {
                inputStream.close();
            }
            if (sc != null) {
                sc.close();
            }
        }
    }

    // given a precursor peptide mass, find the correct isolation window within collection)
    public IsolationWindow FindWindow(double mass, int charge) {
        double mz = Utils.MassToMz(mass, charge);
        int len = windows.size();
        int found = FindWindowIndex(mz, len - 1, 0);
        if (found >= 0) {
            return windows.get(found);
        }
        return null;
    }

    // given mz, return correct isolation window index
    // binary search the window that the m/z could be
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

}

 */
