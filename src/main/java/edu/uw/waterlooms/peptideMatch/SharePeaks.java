package edu.uw.waterlooms.peptideMatch;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class SharePeaks {
    public static ArrayList<PSMResult> ReadPSMFromResultFile(String resultPSMFileName) throws IOException {
        ArrayList<PSMResult> listPSMResult = new ArrayList<>();

        FileReader reader;
        try {
            reader = new FileReader(resultPSMFileName);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(reader)) {
            String line;

            br.readLine();//skip first line
            while ((line = br.readLine()) != null) {
                String[] range = line.split("\\s+");


                PSMResult psm = new PSMResult(
                        range[0],
                        Long.parseLong(range[1]),
                        range[2],
                        range[3],
                        Double.parseDouble(range[4]),
                        Double.parseDouble(range[5]),
                        Double.parseDouble(range[6]),
                        Double.parseDouble(range[7]),
                        Integer.parseInt(range[8]),
                        "",
                        "",
                        Double.parseDouble(range[9]),
                        Double.parseDouble(range[10]),
                        Double.parseDouble(range[11])

                );
                listPSMResult.add(psm);
            }

            br.close();
            reader.close();
        }
        return listPSMResult;

    }
}
