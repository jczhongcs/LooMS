package edu.uw.waterlooms.peptideMatch;


import java.util.HashMap;

/*
Residue	3-
letter
code	1-
letter
code	Mono-
isotopic
mass	Average
mass	Structure
Alanine
C3H5NO	Ala	A	71.037114	71.0779	Ala
Arginine
C6H12N4O	Arg	R	156.101111	156.1857	Arg
Asparagine
C4H6N2O2	Asn	N	114.042927	114.1026	Asn
Aspartic acid
C4H5NO3	Asp	D	115.026943	115.0874	Asp
Asn or Asp	Asx	B
Cysteine
C3H5NOS	Cys	C	103.009185	103.1429	Cys
Glutamic acid
C5H7NO3	Glu	E	129.042593	129.114	Glu
Glutamine
C5H8N2O2	Gln	Q	128.058578	128.1292	Gln
Glu or Gln	Glx	Z
Glycine
C2H3NO	Gly	G	57.021464	57.0513	Gly
Histidine
C6H7N3O	His	H	137.058912	137.1393	His
Isoleucine
C6H11NO	Ile	I	113.084064	113.1576	Ile
Leucine
C6H11NO	Leu	L	113.084064	113.1576	Leu
Lysine
C6H12N2O	Lys	K	128.094963	128.1723	Lys
Methionine
C5H9NOS	Met	M	131.040485	131.1961	Met
Phenylalanine
C9H9NO	Phe	F	147.068414	147.1739	Phe
Proline
C5H7NO	Pro	P	97.052764	97.1152	Pro
Serine
C3H5NO2	Ser	S	87.032028	87.0773	Ser
Threonine
C4H7NO2	Thr	T	101.047679	101.1039	Thr
Selenocysteine
C3H5NOSe	Sec	U	150.95363	150.0379	SeC
Tryptophan
C11H10N2O	Trp	W	186.079313	186.2099	Trp
Tyrosine
C9H9NO2	Tyr	Y	163.06332	163.1733	Tyr
Unknown	Xaa	X
Valine
C5H9NO	Val	V	99.068414	99.1311	Val
 */
public class AminoAcid {
    // may need to double check mass adjustments
//    https://newbsrcmascot.st-andrews.ac.uk/mascot/help/aa_help.html
//    Amino acid reference data
//    The data in this table are for amino acid residues.
//
//    To calculate the mass of a neutral peptide or protein,
//    sum the residue masses plus the masses of the terminating groups
//            (e.g. H at the N-terminus and OH at the C-terminus).
    public static HashMap<Character, Double> mass = new HashMap<Character, Double>();
    public static double ProtonMass = 1.00727647;
    public static double nTerminusAdjM = -89.029920;
    public static double MOxidization = 15.994915;
    static {
        mass.put('A', 71.037114);
        mass.put('R', 156.101111);
        mass.put('N', 114.042927);
        mass.put('D', 115.026943);
        mass.put('C', 160.030649);//103.009185);//160.03065 add fixed ptm 二琉键 57.021464
        mass.put('Q', 128.058578);
        mass.put('E', 129.042593);
        mass.put('G', 57.021464);
        mass.put('H', 137.058912);
        mass.put('I', 113.084064);
        mass.put('L', 113.084064);
        mass.put('K', 128.094963);
        mass.put('M', 131.040485);
        mass.put('F', 147.068414);
        mass.put('P', 97.052764);
        mass.put('S', 87.032028);
        mass.put('T', 101.047679);
        mass.put('W', 186.079313);
        mass.put('Y', 163.063332);
        mass.put('V', 99.068414);
//        mass.put('X', 196.995499);

    }


    public static double getAminoAcidMass(char cur, char next) {
        double res = mass.get(cur);

        return res;
    }
    public static double getAminoAcidMass(char cur) {
        double res = mass.get(cur);

        return res;
    }

}