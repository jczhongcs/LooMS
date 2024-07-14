package edu.uw.waterlooms;

import com.sun.org.apache.xpath.internal.operations.Or;
import edu.uw.waterlooms.peptideMatch.Config;
import me.tongfei.progressbar.ProgressBar;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.*;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.ArrayList;

public class GenerateTargetDecoyFastaFile {

    public static void main(String[] args) throws IOException {


        Path path = FileSystems.getDefault().getPath("").toAbsolutePath();
        String DOCKER_WORKING_DIR = "/data/";
        String filepath = "";


        String fastaInFile = args[0];
        String fastaInFileRaw = fastaInFile.replaceFirst("[.][^.]+$", "");

//        if ((fastaInFile != null && !fastaInFile.isEmpty())) {
//            // TODO: DIA-WEBAPP submits mzXMLInFile as a path/file.mzXML
//            // TODO: Need to strip this and set rawFileName
//            filepath = fastaInFile.substring(0, fastaInFile.lastIndexOf("/"));
//            DOCKER_WORKING_DIR = filepath + "/";
//        }


        //follow mzxmlFile path set the docker_working_path

        System.out.println(DOCKER_WORKING_DIR);

//        java -jar deBruijn-master/deBruijn.jar --input MQ-e1k2DBD-homo_sapiens_targetProtein.fasta --output concat-MQ-e1k2DBD-homo_sapiens_targetProtein_rep1.fasta
        String [] strDecoyfastafiles =new String[3];

        for(int i = 0 ;i <3;i++)
        {
            ArrayList<String> processBuilderCommand = new ArrayList<String>();
            ArrayList<String> nnCommand = new ArrayList<>(processBuilderCommand);
            nnCommand.add("java");
            nnCommand.add("-jar");
            nnCommand.add("/waterlooms/deBruijn.jar");
            nnCommand.add("--input");
            nnCommand.add(fastaInFile);
            nnCommand.add("--output");
            nnCommand.add( "_rep"+i +".fasta");
            strDecoyfastafiles[i] =  "_rep"+i +".fasta";

            ProcessBuilder processBuilder = new ProcessBuilder();

            processBuilder.command(nnCommand);
            System.out.println(nnCommand);
            try {
                Process process = processBuilder.start();
                int exitCode = process.waitFor();
                System.out.println("\ndeBruijn.jar exited with code : " + exitCode);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }


        FileReader reader;

//        FileWriter fileWriterTest = new FileWriter(Config.fastaFolder+"Combine-concat-MQ-e1k2DBD-Mouse2DecoyFiles.fasta");
        FileWriter fileWriterTest = new FileWriter(fastaInFileRaw+"_1Target3DecoyFiles.fasta");
//        FileWriter fileWriterForDiannTest = new FileWriter(fastaInFileRaw+"_1DecoyR1TargetForDiann.fasta");
        BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterTest);
//        BufferedWriter bufferedWriterDiannTest  = new BufferedWriter(fileWriterForDiannTest);
        String []OriginfastaInfile =new String[2];
        OriginfastaInfile[0] = fastaInFile;
        OriginfastaInfile[1] = strDecoyfastafiles[0] ;


        for(String fastFile:OriginfastaInfile) {
            try {
                reader = new FileReader(fastFile);
            } catch (FileNotFoundException noFile) {
                throw new FileNotFoundException();
            }

            try (BufferedReader br = new BufferedReader(reader)) {
                String line;

                while ((line = br.readLine()) != null) {
                    bufferedWriterTest.write(line + "\n");
                }
                br.close();
                reader.close();
            }
        }
        String []Decoyfastafiles =new String[2];
        Decoyfastafiles[0] = strDecoyfastafiles[1] ;
        Decoyfastafiles[1] = strDecoyfastafiles[2] ;
        for(String fastFile:Decoyfastafiles) {
            try {

                reader = new FileReader(fastFile);
            } catch (FileNotFoundException noFile) {
                throw new FileNotFoundException();
            }

            ProgressBar pb = new ProgressBar("Progress", 930680);
            pb.start();
            try (BufferedReader br = new BufferedReader(reader)) {
                String line;
                StringBuilder curComp = new StringBuilder("");
                String lastTitle = "";

                boolean first = true;
                while ((line = br.readLine()) != null) {
                    pb.step();
                    // if line starts with >: new genome
                    if (line.charAt(0) == '>') {
                        if (first) {
                            lastTitle = line;
                            first = false;
                        } else {
                            // finish up old one
                            bufferedWriterTest.write(lastTitle.replace("DeBruijn","DeBruijn_Decoy") + "\n");
                            bufferedWriterTest.write(curComp + "\n");
                            lastTitle = line;
                            curComp.setLength(0);
                        }

                    } else {
                        curComp.append(line.trim());
                    }
                }
                // add last one
                bufferedWriterTest.write(lastTitle.replace("DeBruijn","DeBruijn_Decoy") + "\n");
                bufferedWriterTest.write(curComp + "\n");

                br.close();
                reader.close();
                pb.stop();
            }
        }
        bufferedWriterTest.flush();
        bufferedWriterTest.close();

        //Generate 1decoy +'R'+noFirstCharOf1Target
        //concat-MQ-_rep0.fasta is output by /waterlooms/deBruijn.jar
        readFastaFileCombineDeBru_R_target("concat-MQ-_rep0.fasta",fastaInFileRaw+"_1DecoyR1TargetForDiann.fasta");

    }



    static void readFastaFileCombineDeBru_R_target(String inFileName, String outFileName) throws IOException {
        String fastaInfile = inFileName;
        FileReader reader;

        FileWriter fileWriterTest = new FileWriter(outFileName);
        BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterTest);

        try {
            reader = new FileReader(fastaInfile);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(reader)) {
            String line;
            StringBuilder curComp = new StringBuilder("");
            String lastTitle = "";
            StringBuilder curDecoy = new StringBuilder("");
            boolean isDecoy = false;

            boolean first = true;
            while ((line = br.readLine()) != null) {
                // if line starts with >: new genome
                if(line.charAt(0) == '>') {
                    if (first) {
                        lastTitle=line;
                        first = false;
                    } else {
                        // finish up old one
                        if(!line.contains("DeBruijn"))
                        {
                            bufferedWriterTest.write(lastTitle+"\n");
                            bufferedWriterTest.write(curDecoy+"R"+curComp.substring(1)+"\n");
                            lastTitle=line;


                            curComp.setLength(0);// = "";
                            curDecoy.setLength(0);// = "";
                            isDecoy = false;

                        }else
                        {
                            isDecoy = true;
                        }
                    }

                } else {
                    if(isDecoy)
                        curDecoy.append(line.trim());
                    else
                        curComp.append(line.trim());
                }
            }
            // add last one
            bufferedWriterTest.write(lastTitle+"\n");
            bufferedWriterTest.write(curDecoy+"R"+curComp.substring(1)+"\n");
            bufferedWriterTest.flush();
            bufferedWriterTest.close();
            br.close();
            reader.close();
        }

    }
}
