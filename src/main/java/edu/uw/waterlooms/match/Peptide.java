package edu.uw.waterlooms.match;

import org.apache.commons.lang3.ArrayUtils;

public class Peptide {
    final double PARTIAL_Y_MASS = 18.0105;

    public String id;
    public String composition;
    public double mass;
    public double mz;
    public int charge;
    public double[] y_ions;
    public double[] b_ions;
    public double[] yIonPlusOne;
    public double[] bIonPlusOne;
    public double[] bIonH2OLoss;
    public double[] bIonNH3Loss;
    public double[] yIonH2OLoss;
    public double[] yIonNH3Loss;

    public Peptide(String id, String composition) {
        this.id = id;
        // calculate mass
        this.composition = composition;
        this.mass = generateMass();
        this.charge = 2; // TODO: Refactor as only charge 2 is accepted currently
        this.mz = MassToMz(this.mass, charge);
        this.GenerateIons();
    }


    public static double MassToMz(double mass, int z) {
        return mass/z + AminoAcid.ProtonMass;
    }

    private double generateMass(){
        int len = composition.length();
        double mass = PARTIAL_Y_MASS;
        for (int i = 0; i < len; i ++){
            char residue = composition.charAt(i);
            mass += AminoAcid.getAminoAcidMass(residue);
        }
        return mass;
    }

    public void GenerateIons() {
        int len = composition.length();
        y_ions = new double[len];
        b_ions = new double[len];
        bIonPlusOne = new double[len];
        yIonPlusOne = new double[len];
        bIonH2OLoss = new double[len];
        bIonNH3Loss = new double[len];
        yIonH2OLoss = new double[len];
        yIonNH3Loss = new double[len];
        double partial_y = 18.0105 ;
        double partial_b = 0; // b has only hydrogen mass in daltons

        int idx = 0;
        while(idx < len){
            char residue = composition.charAt(idx);
            partial_b += AminoAcid.getAminoAcidMass(residue);
            b_ions[idx] = MassToMz(partial_b, this.charge);
            partial_y = mass - partial_b;
            y_ions[idx] = MassToMz(partial_y, this.charge);
            idx ++;
        }
        /*
        Iterator residueIterator = ucsdPeptide.iterator();
        while(residueIterator.hasNext())
        {
            edu.ucsd.msjava.msutil.AminoAcid residue = (edu.ucsd.msjava.msutil.AminoAcid) residueIterator.next();
            partial_b += residue.getAccurateMass();
            b_ions[idx] = MassToMz(partial_b, this.charge);
            partial_y = mass - partial_b;
            y_ions[idx] = MassToMz(partial_y, this.charge);
            idx ++;
        }
        */


        ArrayUtils.reverse(y_ions);
        y_ions[0] = 0;
        b_ions[0] = 0;


        // Forloop set call instead of doing a multi stream
        for (int i = 0; i < len; i++){
            yIonH2OLoss[i] = y_ions[i] - 18.0153;
            yIonNH3Loss[i] = y_ions[i] - 17.02654;
            bIonH2OLoss[i] = b_ions[i] - 18.0153;
            bIonNH3Loss[i] = b_ions[i] - 17.02654;
            yIonPlusOne[i] = y_ions[i] + 1.00727647;
            bIonPlusOne[i] = b_ions[i] + 1.00727647;
        }

        // First y ion and last b ion should be 0 so they do NOT induce spurious matches
        y_ions[0] = 0;
        yIonH2OLoss[0] = 0;
        yIonNH3Loss[0] = 0;
        yIonPlusOne[0] = 0;
        b_ions[0] = 0;
        bIonH2OLoss[0] = 0;
        bIonNH3Loss[0] = 0;
        bIonPlusOne[0] = 0;
    }

    public double getMz(){
        return this.mz;
    }
    public double getMass(){ return this.mass; }

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
}
