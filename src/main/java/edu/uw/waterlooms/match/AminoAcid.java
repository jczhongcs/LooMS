package edu.uw.waterlooms.match;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class AminoAcid {
    // may need to double check mass adjustments
    public static HashMap<Character, Double> mass = new HashMap<Character, Double>();
    public static double ProtonMass = 1.0073;
    public static double nTerminusAdjM = -89.029920;
    public static double MOxidization = 15.994915;
    static {
        mass.put('A', 71.03711);
        mass.put('R', 156.10111);
        mass.put('N', 114.04293);
        mass.put('D', 115.02694);
        mass.put('C', 160.03065);
        mass.put('Q', 128.05858);
        mass.put('E', 129.04259);
        mass.put('G', 57.02146);
        mass.put('H', 137.05891);
        mass.put('I', 113.08406);
        mass.put('L', 113.08406);
        mass.put('K', 128.09496);
        mass.put('M', 131.04049);
        mass.put('F', 147.06841);
        mass.put('P', 97.05276);
        mass.put('S', 87.03203);
        mass.put('T', 101.04768);
        mass.put('W', 186.07931);
        mass.put('Y', 163.06333);
        mass.put('V', 99.06841);
    }

//    public static Set<Character> nTerminatorResidue = new HashSet<>();
//    static {
//        nTerminatorResidue.add('A');
//        nTerminatorResidue.add('C');
//        nTerminatorResidue.add('G');
//        nTerminatorResidue.add('P');
//        nTerminatorResidue.add('S');
//        nTerminatorResidue.add('T');
//        nTerminatorResidue.add('V');
//    }

    public static double getAminoAcidMass(char cur) {
        return mass.get(cur);
    }

}