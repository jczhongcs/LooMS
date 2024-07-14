package edu.uw.waterlooms;

import edu.uw.waterlooms.msutil.OpenMzxml;
import edu.uw.waterlooms.peptideMatch.Genome;
import edu.uw.waterlooms.peptideMatch.Peptide;
import edu.uw.waterlooms.peptideMatch.Utils;
import org.apache.commons.lang3.time.StopWatch;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

public class GenerateAutoRTAndPrositInputFile extends GenerateSpectrumLibraryFile{


    /*
    arg[0] fasta file
    arg[1] mzxml file
     */
    public static void main(String[] args) throws IOException {


        Path path = FileSystems.getDefault().getPath("").toAbsolutePath();
        String DOCKER_WORKING_DIR = "/data/";
        String filepath = "";


        String fastaInFile = args[0];
        String fastaInFileRaw = fastaInFile.replaceFirst("[.][^.]+$", "");

        if ((fastaInFile != null && !fastaInFile.isEmpty())) {
            // TODO: DIA-WEBAPP submits mzXMLInFile as a path/file.mzXML
            // TODO: Need to strip this and set rawFileName
            filepath = fastaInFile.substring(0, fastaInFile.lastIndexOf("/"));
            DOCKER_WORKING_DIR = filepath + "/";
        }


        //follow mzxmlFile path set the docker_working_path

        System.out.println(DOCKER_WORKING_DIR);


        StopWatch peptideGenerateStopWatch = new StopWatch();
        peptideGenerateStopWatch.start();
        ArrayList<Genome> genomes = edu.uw.waterlooms.peptideMatch.FastaFile.ReadFile(fastaInFile);

        // get unique set of peptides:
        Collections.reverse(genomes); // put real ones on top

        List<Peptide> pepsOri = genomes.stream().map(x -> x.arrPeps).flatMap(Collection::stream).filter(x -> x.composition.length() >= edu.uw.waterlooms.peptideMatch.Utils.thresholdPeptideSizeMin
                &&  x.composition.length() <= edu.uw.waterlooms.peptideMatch.Utils.thresholdPeptideSizeMax
        ).sorted().collect(Collectors.toCollection(ArrayList::new));

        if (edu.uw.waterlooms.peptideMatch.Utils.bRemovePeptideFromBothTargetAndDecoy) {

            TreeMap<String, List<Peptide>> mapPepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(edu.uw.waterlooms.peptideMatch.Peptide::getComposition, TreeMap::new, Collectors.toList()));

            TreeMap<String, List<edu.uw.waterlooms.peptideMatch.Peptide>> mapILreplacePepstrListPep = pepsOri.stream().sorted()
                    .collect(Collectors.groupingBy(edu.uw.waterlooms.peptideMatch.Peptide::getReplaceILComposition, TreeMap::new, Collectors.toList()));

            System.out.println(mapPepstrListPep.size());
            ArrayList<edu.uw.waterlooms.peptideMatch.Peptide> arrPepsRemove = new ArrayList<>();

            for (String strPep : mapPepstrListPep.keySet()) {
                boolean bDecoyType = false;
                boolean bTargetType = false;

                if (mapPepstrListPep.get(strPep).size() > 1) {
                    for (edu.uw.waterlooms.peptideMatch.Peptide pep : mapPepstrListPep.get(strPep)) {
                        boolean isDecoy = pep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }
                }
                //REMOVE I L MUTATE peptide in both decoy and target
                if (strPep.contains("I") || strPep.contains("L")) {
                    for (edu.uw.waterlooms.peptideMatch.Peptide muSamePep : mapILreplacePepstrListPep.get(strPep.replaceAll("L", "I"))) {
                        boolean isDecoy = muSamePep.id.contains("DeBruijn");
                        bDecoyType = bDecoyType || isDecoy;
                        bTargetType = bTargetType || !isDecoy;
                    }
                }
                if (!(bDecoyType && bTargetType)) {
                    arrPepsRemove.addAll(mapPepstrListPep.get(strPep));
                }
            }
            System.out.println("The number of arrPepsRemove is: " + arrPepsRemove.size());

            pepsOri = arrPepsRemove;

        }
        List<edu.uw.waterlooms.peptideMatch.Peptide> peps = pepsOri.stream().distinct().sorted().collect(Collectors.toCollection(ArrayList::new));


        System.out.println("The number of peps is: " + peps.size());
        peptideGenerateStopWatch.stop();
        System.out.println("Peptide Generate Elapsed Time in Minutes: " + peptideGenerateStopWatch.getTime(TimeUnit.MINUTES) + " ...");


        String psmOutfile =fastaInFile.replaceFirst("[.][^.]+$", "");

        FileWriter fileWriter = new FileWriter(psmOutfile + "_forAutoRT.csv");
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
        bufferedWriter.write("x\ty\n");
        for(edu.uw.waterlooms.peptideMatch.Peptide p:peps)
        {
            bufferedWriter.write(p.composition+"\t0\n");

        }
        bufferedWriter.flush();
        bufferedWriter.close();
        fileWriter.close();



//        OpenMzxml of = new OpenMzxml(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML");
        OpenMzxml of = new OpenMzxml(args[1]);

        //scan and rt
        //1 ms1 with 22 child ms2 window
        //reference mzvalue intensity
        List<Integer> intKeysList = new ArrayList<>(of.num2scan2.keySet());
        List<Integer> intKeysListMS1 = new ArrayList<>(of.num2scan.keySet());
        double dmzMin  = -1;
        double dmzMax  = -1;
        for(Integer ims1:intKeysListMS1)
        {
            if(dmzMin<0)
            {
                dmzMin = of.num2scan.get(ims1).getScanMzWindowLower();
                dmzMax = of.num2scan.get(ims1).getScanMzWindowUpper();

            }else
            {
                if(dmzMin>of.num2scan.get(ims1).getScanMzWindowLower())
                {
                    dmzMin = of.num2scan.get(ims1).getScanMzWindowLower();

                }
                if(dmzMax<of.num2scan.get(ims1).getScanMzWindowLower())
                {
                    dmzMax = of.num2scan.get(ims1).getScanMzWindowLower();

                }
            }

        }

        double dActiveEnergyLo = of.num2scan2.get(intKeysList.get(0)).getPrecursor().getActivationInfo().getActivationEnergyLo();


        FileWriter forPrositfileWriter = new FileWriter(psmOutfile + "_forProsit.csv");
        BufferedWriter forPorsitbufferedWriter = new BufferedWriter(forPrositfileWriter);
        forPorsitbufferedWriter.write("x\ty\n");
        for(edu.uw.waterlooms.peptideMatch.Peptide p:peps)
        {
            if(p.composition.length() > 30 || p.composition.length()<7) continue;
            for(int i =1 ;i<7;i++)
            {
                double dmz= Utils.MassToMz(p.mass,i);
                if(dmz<dmzMin || dmz>dmzMax) continue;

                forPorsitbufferedWriter.write(p.composition+"\t"+dActiveEnergyLo+"\t"+i+"\n");

            }


        }
        forPorsitbufferedWriter.flush();
        forPorsitbufferedWriter.close();
        forPrositfileWriter.close();

    }

}
