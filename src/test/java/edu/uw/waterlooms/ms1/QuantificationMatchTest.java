package edu.uw.waterlooms.ms1;

import javolution.io.Struct;
import org.junit.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

enum COMPARE_TYPE  {
    MS1_TRAIL,
    FEATURE
}
public class QuantificationMatchTest {
    private COMPARE_TYPE compare_type;
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
    public static void main(String[] args) throws IOException {
        String current_path = System.getProperty("user.dir");
        String fname1 = "R01_ms1_trails.tsv";
        String fname2 = "R02_ms1_trails.tsv";
        String fpath1 = current_path + "/data/" + fname1;
        String fpath2 = current_path + "/data/" + fname2;
        double ppm = 5e-6;
        double rt_tol = 0.5;

        QuantificationMatchTest quantificationMatchTest = new QuantificationMatchTest();
        quantificationMatchTest.compare_type = COMPARE_TYPE.MS1_TRAIL;
//        quantificationMatchTest.compare_type = COMPARE_TYPE.FEATURE;

        if (quantificationMatchTest.compare_type == COMPARE_TYPE.MS1_TRAIL) {
            quantificationMatchTest.mz = 0;
            quantificationMatchTest.rt = 1;
            quantificationMatchTest.peaks_sum = 2;
            quantificationMatchTest.peaks_area = 3;
        }
        quantificationMatchTest.intersectList = new ArrayList<>();
        quantificationMatchTest.exceptionList = new ArrayList<>();
//        quantificationMatchTest.get_union();

        InOutFile read_isoList1 = new InOutFile();
        InOutFile read_isoList2 = new InOutFile();
        read_isoList1.inputFile(fpath1);
        read_isoList2.inputFile(fpath2);
        quantificationMatchTest.init(read_isoList1.list, read_isoList2.list);
        quantificationMatchTest.intersect(ppm,rt_tol,quantificationMatchTest.isoList1,quantificationMatchTest.isoList2);


        InOutFile ms1_write_intersection = new InOutFile(quantificationMatchTest.intersectList);
        ms1_write_intersection.writeFile(current_path + "/data/" + "R01_peak_ms1_trails_R02.tsv",
                "mz_R01\trt_R01\tpeaks_sum_R01\tmz_R01\trt_R01\tpeaks_sum_R01");
        InOutFile ms1_write_exception = new InOutFile(quantificationMatchTest.exceptionList);
        ms1_write_exception.writeFile(current_path + "/data/" + "R01_exception_peak_ms1_trails_R02.tsv",
                "mz\trt\tpeaks_sum_R01\tpeaks_area_R01\tquality_R01");
//        ms1_write.writeFile(current_path + "/data/" + "R01_peaks_sne",
//                "index\tprecursor_mz\trt\tassnum\tpepseq\tpeaks_sum\tpeaks_area");
    }

    private void init(List<Double[]> list1, List<Double[]> list2) {
        isoList1 = list1;
        isoList2 = list2;
        isoList1.sort(Comparator.comparing(l -> l[rt]));
        isoList2.sort(Comparator.comparing(l -> l[rt]));
        if (compare_type == COMPARE_TYPE.FEATURE) {
            isoList1.sort(Comparator.comparing(l -> l[z]));
            isoList2.sort(Comparator.comparing(l -> l[z]));
        }
        isoList1.sort(Comparator.comparing(l -> l[mz]));
        isoList2.sort(Comparator.comparing(l -> l[mz]));
    }

    public void intersect(double ppm, double rt_tol, List<Double[]> isoList1, List<Double[]> isoList2){
        Boolean[] occupiedList2 = new Boolean[isoList2.size()];
        for (int i = 0; i < occupiedList2.length; i++) {
            occupiedList2[i] = false;
        }

        int j = 0;
        for (Double[] iso1: isoList1) {
            double mz1 = iso1[mz];
            double rt1 = iso1[rt];
            boolean if_found = false;
            ArrayList<Integer> allMatch = new ArrayList<>();

            while (j < isoList2.size() && isoList2.get(j)[mz] > mz1 * (1 - ppm)) {
                if (j == 0) break;
                j--;
            }
            for (; j < isoList2.size(); j++) {
                Double[] iso2 = isoList2.get(j);
                double rt2 = iso2[rt];
                double mz2 = iso2[mz];
                if (mz2 >= mz1 * (1 - ppm) && mz2 <= mz1 * (1 + ppm)
                        && rt2 >= rt1 - rt_tol && rt2 <= rt1 + rt_tol
                        && !occupiedList2[j]) {
                    if ((this.compare_type == COMPARE_TYPE.FEATURE && iso1[z].intValue() == iso2[z].intValue())
                    || this.compare_type == COMPARE_TYPE.MS1_TRAIL) {
                        allMatch.add(j);
                        if_found = true;
                    }
                }
                if (mz2 > mz1 * (1 + ppm)) {
                    break;
                }
            }
            if (if_found) {
                // ppm difference the smallest
                int pos = 0;
                double minDifference = Math.abs((mz1 - isoList2.get(allMatch.get(pos))[mz])/mz1);
                for (int d = 1; d < allMatch.size(); d++) {
                    if (minDifference > Math.abs((mz1 - isoList2.get(allMatch.get(d))[mz])/mz1) ) {
                        pos = d;
                        minDifference = Math.abs((mz1 - isoList2.get(allMatch.get(d))[mz])/mz1);
                    }
                }
                Double[] iso2 = isoList2.get(allMatch.get(pos));
                if (this.compare_type == COMPARE_TYPE.FEATURE) {
                    Double[] newIso = new Double[8];
                    newIso[0] = iso1[mz];
                    newIso[1] = iso1[rt];
                    if (iso1.length > 15 && iso2.length > 15) {
                        newIso[2] = iso1[peaks_sum];
                        newIso[3] = iso1[peaks_area];
                        newIso[5] = iso2[peaks_sum];
                        newIso[6] = iso2[peaks_area];
                        newIso[7] = iso2[quality_score];
                    }
                    intersectList.add(newIso);
                } else {
                    Double[] newIso = new Double[9];
                    newIso[0] = iso1[mz];
                    newIso[1] = iso1[rt];
                    newIso[2] = iso1[peaks_sum];
                    newIso[3] = iso1[peaks_area];
                    newIso[4] = iso2[mz];
                    newIso[5] = iso2[rt];
                    newIso[6] = iso2[peaks_sum];
                    newIso[7] = iso2[peaks_area];
                    newIso[8] = iso1[rt] - iso2[rt];
                    intersectList.add(newIso);
                }
                occupiedList2[allMatch.get(pos)] = true;
            } else { // not found
                if (this.compare_type == COMPARE_TYPE.FEATURE) {
                    Double[] newIso = new Double[5];
                    newIso[0] = iso1[mz];
                    newIso[1] = iso1[rt];
                    if (iso1.length > 15 && iso1.length > 15) {
                        newIso[2] = iso1[peaks_sum];
                        newIso[3] = iso1[peaks_area];
                        newIso[4] = iso1[quality_score];
                    }
                    exceptionList.add(newIso);
                } else {
                    Double[] newIso = new Double[4];
                    newIso[0] = iso1[mz];
                    newIso[1] = iso1[rt];
                    newIso[2] = iso1[peaks_sum];
                    newIso[3] = iso1[peaks_area];
                    exceptionList.add(newIso);
                }
            }
        }
        double total = isoList1.size();
        double percentage_intersectList = intersectList.size()/total * 100;
        double percentage_exceptionList = exceptionList.size()/total * 100;
        System.out.println("ppm: " + ppm + ", rt: " + rt_tol);
        System.out.println(percentage_intersectList + "% (" + intersectList.size() + "/" + total + ")");
        System.out.println(percentage_exceptionList + "% (" + exceptionList.size() + "/" + total + ")");
    }

    public void get_union() throws IOException {
        String current_path = System.getProperty("user.dir");
        String fname1 = "R03_exception_peak_R01.tsv";
        String fname2 = "R03_exception_peak_R02.tsv";
        String output_file = "R03_exception_peak.tsv";
        String fpath1 = current_path + "/data/" + fname1;
        String fpath2 = current_path + "/data/" + fname2;

        InOutFile read_isoList1 = new InOutFile();
        InOutFile read_isoList2 = new InOutFile();
        read_isoList1.inputFile(fpath1);
        read_isoList2.inputFile(fpath2);
        List<Double[]> iso1;
        List<Double[]> iso2;
        iso1 = read_isoList1.list;
        iso2 = read_isoList2.list;
        List<Double[]> union;

        iso1.sort(Comparator.comparing(l -> l[1]));
        iso2.sort(Comparator.comparing(l -> l[1]));
        iso1.sort(Comparator.comparing(l -> l[0]));
        iso2.sort(Comparator.comparing(l -> l[0]));
        union = new ArrayList<>();

        int common = 0;
        int num_first = 0;
        int num_second = 0;
        int i = 0;
        int j = 0;
        while (i < iso1.size()) {
            if (j == iso2.size()) {
                break;
            }
            Double[] iso_first = iso1.get(i);
            double mz1 = iso_first[0];
            double rt1 = iso_first[1];
            Double[] iso_second = iso2.get(j);
            double mz2 = iso_second[0];
            double rt2 = iso_second[1];
            if (mz2 == mz1  && rt1 == rt2) {
                Double[] newIso = new Double[5];newIso[0] = iso_first[0];
                    newIso[1] = iso_first[1];
                    newIso[2] = iso_first[2];
                    newIso[3] = iso_first[3];
                    newIso[4] = iso_first[4];
                    union.add(newIso);
                    common++;
                    i++;
                    j++;
                } else if ((mz2 == mz1  && rt1 > rt2) || mz2 < mz1) {
                    Double[] newIso = new Double[5];
                    newIso[0] = iso_second[0];
                    newIso[1] = iso_second[1];
                    newIso[2] = iso_second[2];
                    newIso[3] = iso_second[3];
                    newIso[4] = iso_second[4];
                    union.add(newIso);
                    num_second++;
                    j++;
                } else if ((mz2 == mz1  && rt1 < rt2) || mz2 > mz1) {
                    Double[] newIso = new Double[5];
                    newIso[0] = iso_first[0];
                    newIso[1] = iso_first[1];
                    newIso[2] = iso_first[2];
                    newIso[3] = iso_first[3];
                    newIso[4] = iso_first[4];
                    union.add(newIso);
                    num_first++;
                    i++;
                } else {
                    System.out.println("This message should not appear.");
                }

                if (j == iso2.size()) {
                    Double[] newIso = new Double[5];
                    newIso[0] = iso_first[0];
                    newIso[1] = iso_first[1];
                    newIso[2] = iso_first[2];
                    newIso[3] = iso_first[3];
                    newIso[4] = iso_first[4];
                    union.add(newIso);
                    num_first++;
                    i++;
                }
        }
        if (i == iso2.size()) {
            while (j <= iso2.size()) {
                Double[] newIso = new Double[5];
                Double[] iso_second = iso2.get(j);
                newIso[0] = iso_second[0];
                newIso[1] = iso_second[1];
                newIso[2] = iso_second[2];
                newIso[3] = iso_second[3];
                newIso[4] = iso_second[4];
                union.add(newIso);
                num_second++;
                j++;
            }
        }
        System.out.println("Common: " + common);
        System.out.println("Only in first: " + num_first);
        System.out.println("Only in second: " + num_second);
        InOutFile write_union = new InOutFile(union);
        write_union.writeFile(current_path + "/data/" + output_file,
                "mz\trt\tpeaks_sum_R01\tpeaks_area_R01\tquality_R01");
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
