/*
package edu.uw.waterlooms.database;


import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.match.IsotopeTrail;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

public class IsolationWindowSerializationService {
    private String xicInfile;

    public IsolationWindowSerializationService(String xicInfile){
        this.xicInfile = xicInfile;
    }

    public void serializeXIC(String outDir) throws IOException {

        FileInputStream inputStream = null;
        Scanner sc = null;
        inputStream = new FileInputStream(this.xicInfile);
        sc = new Scanner(inputStream, "UTF-8");

        int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail
        int windowCounter = 0;

        IsolationWindow cur;
        IsotopeTrail curTrail = new IsotopeTrail();
        ArrayList<IsolationWindow> windows =  new ArrayList();

        while(sc.hasNextLine()){
            String line = sc.nextLine();
            if (line.isEmpty()) { continue; }
            if (line.contains("END")){
                // Sort Trails Here
                Collections.sort(cur.trails);

                // Serialize the Isolation window to n_window.ser
                FileOutputStream xicOut = new FileOutputStream(outDir + windowCounter + "_window.ser");
                ObjectOutputStream out = new ObjectOutputStream(xicOut);
                out.writeObject(cur);
                out.close();
                xicOut.close();

                // Add the value(s) to a lookup table here
                //windows.add(cur); Do NOT append to a list
                continue;
            }
            if (line.contains("START")) {
                // For Progress Tracking
                System.out.println(windowCounter + "/" + "22");
                windowCounter += 1;
                cur = new IsolationWindow();
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
                curTrail.mz = Double.parseDouble(mzLine[0]);
                curTrail.maxIntensityRt = Double.parseDouble(mzLine[1]);
                lineType = 1;
                continue;
            }
            if (lineType == 1) {
                String[] linelst = line.split("\\s+");
                curTrail.ParseIntensities(linelst);
                lineType = 2;
                continue;
            }
            if (lineType == 2) {
                String[] linelst = line.split("\\s+");
                curTrail.ParseRts(linelst);
                // add to list
                cur.PutRtMax(curTrail.rts, curTrail.intensities);
                cur.trails.add(curTrail);
                curTrail = new IsotopeTrail();
                lineType = 0;
                continue;
            }
        }
    }

}
*/