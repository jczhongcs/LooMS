package edu.uw.waterlooms.match;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class FastaFile {
    // taskes in a fasta file and parse into genomes
    public static ArrayList<Genome> ReadFile(String fastaInfile) throws IOException {
        ArrayList<Genome> res =  new ArrayList<>();
        // String dataFile = file;

        FileReader reader;
        try {
            reader = new FileReader(fastaInfile);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(reader)) {
            String line;
            String curComp = "";
            String curId = "";
            String curDesc = "";
            boolean first = true;
            while ((line = br.readLine()) != null) {
                // if line starts with >: new genome
                if(line.charAt(0) == '>') {
                    if (first) {
                        first = false;
                    } else {
                        // finish up old one
                        res.add(new Genome(curId,curComp,curDesc));
                    }

                    // start new one
                    String[] lineLst = line.split("\\|");
                    curDesc = lineLst[lineLst.length -1];
                    curId = String.join("" , Arrays.copyOfRange(lineLst, 0, lineLst.length - 1));
                    curComp = "";
                } else {
                    curComp+=line.trim();
                }
            }
            // add last one
            res.add(new Genome(curId,curComp,curDesc));
            br.close();
            reader.close();
        }
        return res;
    }
}
