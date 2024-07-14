package edu.uw.waterlooms.ms2;

//import javolution.io.Struct;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class MS2FeatureQuantificationTest {
    public int id = 0;
    public int mz = 1;
    public int rt = 2;
    private int z = 3;
    private int quality_score = 14;
    private int peaks_sum = 15;
    private int peaks_area = 16;
    List<Double[]> isoList1;
    List<Double[]> isoList2;
    List<Double[]> intersectList;
    List<Double[]> exceptionList;

//    @Test
    public void getPeptideSequenceWithMZAndRT() throws IOException {
        String current_path = System.getProperty("user.dir");
        String name = "spectranaut_results_number_only.tsv";
        String spectranautFilePath = current_path + "/data/" + name;
        InOutFile spectranautFile = new InOutFile();
        spectranautFile.inputFile(spectranautFilePath);
        double mzTolPPM = 5e-6;
        double rtTol = 0.2;
        List<Double[]> results = new ArrayList<>();

        for (int i = 0; i < spectranautFile.list.size();) {
            Double[] l1 = spectranautFile.list.get(i);
            double mz1 = l1[0];
            double rt1 = l1[1];
            double z1 = l1[2];


            if (i+1<spectranautFile.list.size() && i+2<spectranautFile.list.size()) {
                Double[] l2 = spectranautFile.list.get(i+1);
                Double[] l3 = spectranautFile.list.get(i+2);
                double mz2 = l2[0];
                double rt2 = l2[1];
                double z2 = l2[2];
                double mz3 = l3[0];
                double rt3 = l3[1];
                double z3 = l3[2];

                if (mz1 <= mz2 * (1 + mzTolPPM) && mz1 >= mz2 * (1 - mzTolPPM)
                        && rt1 <= (rt2 + rtTol) && rt1 >= (rt2 - rtTol)
                        && z1 == z2) {
                    results.add(l1);
                    i++;
                    if (mz1 <= mz3 * (1 + mzTolPPM) && mz1 >= mz3 * (1 - mzTolPPM)
                            && rt1 <= (rt3 + rtTol) && rt1 >= (rt3 - rtTol)
                            && z1 == z3) {
                        i++;
                    }
                }
            }
            i++;
        }
        spectranautFile.list = results;
        String outputPath = current_path + "/data/single_peptide.tsv";
        spectranautFile.writeFile(outputPath, "FG.PrecMz\tEG.ApexRT\tFG.Charge\tFG.MS1Quantity\tFG.Quantity");
    }


}

class InOutFile {
    List<Double[]> list;
    InOutFile() {}
    InOutFile(List<Double[]> list) {
        this.list = list;
    }
    public void inputFile(String filename) throws IOException {
        List<Double[]> result = new ArrayList<>();
        FileReader fileReader = new FileReader(filename);
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        String line;
        boolean IsfirstLine = true;
        while((line = bufferedReader.readLine()) != null) {
            // process first line
            if (IsfirstLine == true) {
                IsfirstLine = false;
            }
            else {
                //process line
                String[] words = line.split("\\t");
                Double[] nums = new Double[words.length];
                for (int i = 0; i < words.length; i++) {
                    nums[i] = Double.valueOf(words[i]);
                }
                result.add(nums);
            }
        }
        this.list = result;
    }

    // simply write the list
    public void writeFile(String filename, String header) throws IOException {
        File file = new File(filename);
        try {
            // create FileWriter object with file as parameter
            FileWriter outputfile = new FileWriter(file);
            // PrintWriter
            PrintWriter printWriter = new PrintWriter(outputfile);
            // add header
            printWriter.print(header + '\n');
            int line_size = list.get(0).length;
            for (int i = 0; i < list.size(); i++) {
                String data = list.get(i)[0].toString();
                for (int j = 1; j < line_size; j++){
                    data += '\t' + list.get(i)[j].toString();
                }
                data += '\n';
                printWriter.print(data);
            }
            printWriter.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
}

