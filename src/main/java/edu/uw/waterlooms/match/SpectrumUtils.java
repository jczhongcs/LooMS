package edu.uw.waterlooms.match;
import edu.uw.waterlooms.entity.IsolationWindow;
import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.Arrays;

public class SpectrumUtils {
    public static void main(String[] args) throws Throwable{
        String [] files = {"window_1300.mgf"};
//        buildSpecCsvForFiles(files);
        buildSpecMgfFromCsv(Config.userHome + Config.spectrumFolder+"sample.mgf", Config.userHome + Config.spectrumFolder+"spectrum.csv");
    }

    // builds csv file from mgf spectrum
    /*
    public static void buildSpecCsvForFiles(String [] files) throws Throwable{
//        String [] files = {"all_trials_ms2_4.mgf"};
        File outputFile= new File(Config.userHome + Config.spectrumFolder+"spectrum_rank.csv");
        FileWriter fw;
        StringBuilder sb = new StringBuilder("\n");

        for (String file:files) {
            IsolationWindowCollection spec = new IsolationWindowCollection(file);
            // create a spectrum csv
            for (IsolationWindow w : spec.windows) {
                sb.append(w.windowToCsv());
            }
        }

        sb.append("\n");
        if (outputFile.exists())
        {
            fw = new FileWriter(outputFile,true); //if file exists append to file. Works fine.
        }
        else
        {
            outputFile.createNewFile();
            fw = new FileWriter(outputFile);
        }
        fw.write(sb.toString());
    }\

     */

    // given mgf and csv with half ranks calculated, add half rank of each peak for every trail after retention times.
    // currently, the code assumes peaks in mgf and in csv are in the same order and appears only once
    // (they must match exactly);
    // this program does not work fully yet as currently, mgf might have duplicated peaks
    // there are caveats with floating point number associated with file reading and rounding; floating point number
    // manipulation done in this function is currently correct, so please be careful when changing anything!
    public static void buildSpecMgfFromCsv(String mgf, String csv) throws IOException {
        DecimalFormat massDf = new DecimalFormat("#.#####");
        massDf.setRoundingMode(RoundingMode.HALF_UP);

        DecimalFormat rtDf = new DecimalFormat("#.############");
        rtDf.setRoundingMode(RoundingMode.HALF_UP);

        File outputFile= new File(Config.userHome + Config.spectrumFolder+"sample_rank.mgf");
        FileWriter fw;
        FileReader mgfReader;
        FileReader csvReader;
        BufferedReader mgfLine;
        BufferedReader csvLine;
        StringBuilder sb = new StringBuilder();
        String csvCur;
        try {
            mgfReader = new FileReader(mgf);
            csvReader = new FileReader(csv);
            mgfLine = new BufferedReader(mgfReader);
            csvLine = new BufferedReader(csvReader);
            csvCur = csvLine.readLine();
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        String line;
        int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail
        String mzLow = "";
        String curMass = "";
        String[] range;
        while ((line = mgfLine.readLine()) != null) {
            sb.append(line + "\n");
            // if line empty, continue
            if (line.isEmpty()) { continue; }

            // if see end, add cur to list, next
            if (line.contains("END")) {
                continue;
            }

            if (line.contains("START")) {
                line = mgfLine.readLine();
                sb.append(line+ "\n");
                range = line.split("\\s+");
                mzLow = range[0];
                lineType = 1;
                continue;
            }

            if (lineType == 1) {
                lineType = 0;
                range = line.split("\\s+");
                curMass = massDf.format(Double.parseDouble(range[0].trim()));
            }

            sb.append(mgfLine.readLine()+ "\n");

            line = mgfLine.readLine();
            range = line.split("\\s+");
            sb.append(line+ "\n");

            lineType = 1;
            // add the line

            String ranks = "";
            String[] csvCurLst = csvCur.split(",");
            boolean mzlowCond = csvCurLst[1].equals(mzLow.trim());
            boolean massCond = massDf.format(Double.parseDouble(csvCurLst[2])).equals(curMass);
            boolean rtCond = Arrays.asList(range).contains(rtDf.format(Double.parseDouble(csvCurLst[3])));

            boolean debug = false;
            // makes sure window, mz and retention time all match to identify the same peak
            while(mzlowCond && massCond && rtCond) {
                debug = true;
                ranks += csvCurLst[csvCurLst.length-1];
                ranks += " ";
                csvCur = csvLine.readLine();
//                if(csvCur==null) { break; }
                csvCurLst = csvCur.split(",");
                mzlowCond = csvCurLst[1].equals(mzLow.trim());
                massCond = massDf.format(Double.parseDouble(csvCurLst[2])).equals(curMass.trim());
                rtCond = Arrays.asList(range).contains(rtDf.format(Double.parseDouble(csvCurLst[3])));
            }
            if (!debug) {
                int a = 0;
            }
            ranks+="\n";
            sb.append(ranks);
        }
        mgfReader.close();
        csvReader.close();
        mgfLine.close();
        csvLine.close();

        if (outputFile.exists())
        {
            fw = new FileWriter(outputFile,false); //if file exists append to file. Works fine.
        }
        else
        {
            outputFile.createNewFile();
            fw = new FileWriter(outputFile);
        }
        fw.write(sb.toString());

    }


}
