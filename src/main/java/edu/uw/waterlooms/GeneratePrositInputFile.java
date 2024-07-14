package edu.uw.waterlooms;

import edu.uw.waterlooms.msutil.OpenMzxml;
import edu.uw.waterlooms.peptideMatch.Config;

import java.util.ArrayList;
import java.util.List;

public class GeneratePrositInputFile extends GenerateSpectrumLibraryFile{
    public static void main(String args[])
    {
        OpenMzxml of = new OpenMzxml(Config.spectrumFolder + "Fig1_MP-DIA-120min120kMS1-22W30k-8dppp_MHRM_R01.mzXML");

        //scan and rt
        //1 ms1 with 22 child ms2 window
        //reference mzvalue intensity
        List<Integer> intKeysList = new ArrayList<>(of.num2scan2.keySet());

        double dActiveEnergyLo = of.num2scan2.get(intKeysList.get(0)).getPrecursor().getActivationInfo().getActivationEnergyLo();


    }
}
