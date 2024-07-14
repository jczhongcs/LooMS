package edu.uw.waterlooms.peptideMatch;


public class Peptide implements Comparable<Peptide>{
    public String id;

    public String getComposition() {
        return composition;
    }
    public String getReplaceILComposition() {
        return composition.replaceAll("L","I");
    }
    public String composition;
    public double mass;
    // y/b ions are populated during DB match, since if peptide is not matched, there's no need to generate them
    public double[] y_ions; // charge 0 mass
    public double[] b_ions;

    public String chemicalComponent;//chemical component
    public double[] ad_bIntensity;//b ions intensity
    public double[] ad_yIntensity;//b ions intensity

    public int iMissCleavage;

    public Peptide(String id, String substring, double v) {
        this(id,substring);
        setdMutationRate(v);

    }

    public void setdMutationRate(double dMutationRate) {
        this.dMutationRate = dMutationRate;
    }

    public double dMutationRate;//相对正常切的变异系数


    public double dPeakAreaSum;//peakArea sum


//    public static double MassToMz(double mass, int z) {
//        return mass/z + AminoAcid.ProtonMass;
//    }
    public Peptide(String id, String composition) {
        this.id = id;
        // calculate mass
        double sum = 18.0105647D;
        int len = composition.length();
        for (int i = 0; i < len; ++i) {
//            char next = i == len -1 ? '\0' : composition.charAt(i+1);
            sum += AminoAcid.getAminoAcidMass(composition.charAt(i));
        }

        this.composition = composition;
        this.mass = sum;
//        GenerateIons(); // generate y/b ions here

    }

    public void GenerateIons() {
        if(b_ions==null) {
            int len = composition.length();
//            int len = composition.length() - 1;
//            double partial_y = 18.0105647D + 1.00728D;
//            double partial_b = 1.00728D;
            double partial_y = 18.0105647D + AminoAcid.ProtonMass;
            double partial_b = AminoAcid.ProtonMass;
            y_ions = new double[len];
            b_ions = new double[len];
            for (int i = 0; i < len; ++i) {
                char next = i == len - 1 ? '\0' : composition.charAt(i + 1);
                partial_y += AminoAcid.getAminoAcidMass(composition.charAt(len - i -1), next);//AminoAcid.mass.get(composition.charAt(len - i));
                partial_b += AminoAcid.getAminoAcidMass(composition.charAt(i), next);//AminoAcid.mass.get(composition.charAt(i));
                b_ions[i] = partial_b; //MassToMz(partial_b, charge);
                y_ions[len - i - 1] = partial_y; //MassToMz(partial_y, charge);
            }
        }
    }

    @Override
    public boolean equals(Object o) {
        if (o instanceof Peptide) {
            return this.composition.compareTo(((Peptide) o).composition) == 0;
        }
        return false;
    }

    @Override
    public int hashCode() {
        return composition.hashCode();
    }

    @Override
    public int compareTo(Peptide o) {
        return Double.compare(mass,o.mass);
    }

    public static double CalculateMassByPeptideStr(String strPep)
    {
        double mass = 0.0;
        double sum = 18.0105647D;
        int len = strPep.length();
        for (int i = 0; i < len; ++i) {
//            char next = i == len -1 ? '\0' : composition.charAt(i+1);
            sum += AminoAcid.getAminoAcidMass(strPep.charAt(i));
        }

        mass = sum;
        return mass;
    }
}
