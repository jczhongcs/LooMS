package edu.uw.waterlooms.peptideMatch;


import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

public class FastaFile {
    // taskes in a fasta file and parse into genomes
    static Map<Character,Integer> mapAnomialAcidFrequence= new HashMap<>();
    static Map<Character,Double> mapAnomialAcidFrequenceRate= new HashMap<>();

    static Integer iCount=0;

    public static ArrayList<Genome> ReadFile(String fastaInfile) throws IOException {
        ArrayList<Genome> res =  new ArrayList<>();
        // String dataFile = file;

        FileReader reader;
        try {
            System.out.println(fastaInfile);
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
                        addResGenome(res, curComp, curId, curDesc);


                        //    System.out.println(res.size());
//                        res.add(new Genome(curId,curComp,curDesc));
                    }

                    // start new one
                    String[] lineLst = line.split("\\|");
                    curDesc = lineLst[lineLst.length -1];
                    if( lineLst.length>1)
                    {
                         curId = String.join("" , Arrays.copyOfRange(lineLst, 0, lineLst.length - 1));
                    }
                    else
                    {
                         curId = String.join("" , Arrays.copyOfRange(lineLst, 0, lineLst.length ));
                    }
//                    curId = String.join("" , Arrays.copyOfRange(lineLst, 0, lineLst.length - 1));
                    curComp = "";
                } else {
                    curComp+=line.trim();
                }
            }
            // add last one
            addResGenome(res, curComp, curId, curDesc);
//            res.add(new Genome(curId,curComp.substring(1),curDesc));
            br.close();
            reader.close();
        }

        tranFrequenceToRate();
        return res;
    }

    private static void addResGenome(ArrayList<Genome> res, String curComp, String curId, String curDesc) {
        if(Utils.FastaReadEnumMode == Utils.FastaReadEnum.TrageMode)
        {
            if(!curId.contains("DeBruijn"))
            {
                res.add(new Genome(curId, curComp.substring(1), curDesc));
            }
        }else if(Utils.FastaReadEnumMode == Utils.FastaReadEnum.DecoyMode)
        {
            if(curId.contains("DeBruijn"))
            {
                res.add(new Genome(curId, curComp.substring(1), curDesc));
            }
        }else{
            res.add(new Genome(curId, curComp.substring(1), curDesc));
        }

        putAnomialFrequence(mapAnomialAcidFrequence,curComp.substring(1),false);

    }

    private static void putAnomialFrequence(Map<Character, Integer> map, String s,boolean bTarget) {
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            Integer val = map.get(c);
            if (val != null) {
                map.put(c, new Integer(val + 1));
            }
            else {
                map.put(c, 1);
            }
            iCount++;

        }
    }

    public static double getFrequenceRateMultiple(String s)
    {
        double dResult = 1.0;
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            dResult *= mapAnomialAcidFrequenceRate.get(c);
        }
        return Math.log(dResult);
    }
    public static void tranFrequenceToRate()
    {
        mapAnomialAcidFrequence.forEach((x,y)->
        {
            mapAnomialAcidFrequenceRate.put(x,(y/(double)iCount));
        });
    }

    public static void outputAnomialAcidFrequency(String fileName)
    {
        FileWriter fileWriter = null;

        try {
            fileWriter =  new FileWriter(fileName);
        } catch (IOException e) {
            e.printStackTrace();
        }
        BufferedWriter bufferedWriter  = new BufferedWriter(fileWriter );

        mapAnomialAcidFrequence.forEach((x,y)->
        {
            try {
                bufferedWriter.write(x+"\t"+(y/(double)iCount)+"\n");
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

        });


        try {
            bufferedWriter.flush();
            bufferedWriter.close();
            fileWriter.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }


    }

    public static ArrayList<Genome> ReadFileFormPeptide(String fastaInfile)  throws IOException{

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
            Genome ge = new Genome(curId,"",curDesc);
            while ((line = br.readLine()) != null) {
                // if line starts with >: new genome
                if(line.charAt(0) == '>') {


                    // start new one
                    String[] lineLst = line.split("\\|");
                    curDesc = lineLst[lineLst.length -1];
                    if( lineLst.length>1)
                    {
                        curId = String.join("" , Arrays.copyOfRange(lineLst, 0, lineLst.length - 1));
                    }
                    else
                    {
                        curId = String.join("" , Arrays.copyOfRange(lineLst, 0, lineLst.length ));
                    }
//                    curId = String.join("" , Arrays.copyOfRange(lineLst, 0, lineLst.length - 1));
                    curComp = "";
                } else {
                    Peptide pe = new Peptide(curId, line.trim());
                    ge.arrPeps.add(pe);
//                    curComp+=line.trim();
                }
            }
            // add last one
            res.add(ge);
            br.close();
            reader.close();
        }
        return res;
    }

    public static ArrayList<Genome> ReadFormPSMList(ArrayList<PSMResult> arrPsmResults) {
        ArrayList<Genome> res =  new ArrayList<>();
        Genome ge = new Genome("","","");

        for (PSMResult psm:arrPsmResults
             ) {
            Peptide pe = new Peptide(psm.proteinname,psm.strpeptide);
            ge.arrPeps.add(pe);
            res.add(ge);

        }
        return res;

    }
}
