package edu.uw.waterlooms.peptideMatch;

import java.util.Objects;

public class PSMResult implements  Comparable<PSMResult> {
    public int sharePeaksAll;
    String strpeptide;
    Long msonefeature;
    String proteinname;
    String decoyTarget;
    double pepmass;
    double msonemass;
    double msonetime;
    double msoneMZ;
    int msoneZ;
//    bionMatchCount_C1\tyionMatchCount_C1\tbionMatchCount_C2\tyionMatchCount_C2\tibMatchNormalAndNH3loss" +
//    "\tiyAllMatchedWithloss\t" +
    String BIon;
    String YIon;

    int iNo;
    double maxShareMS1MZ ;
    double maxShareMS1Time;
    Long ims1id;


    public double getTargetscore() {
        return targetscore;
    }

    //    \tdecoyLabel\ttargetlabel\tdecoyscore\t
    double targetscore;
    double countDP;
    double fdr;
    int sharepeaksMax;

    String strSharePeaksMaxInfo;

    String strSharePeaksAllInfo;

    public PSMResult(String strpeptide, Long msonefeature, String proteinname, String decoyTarget,
                     double pepmass, double msonemass, double msonetime, double msoneMZ, int msoneZ,
                     String BIon, String YIon, double targetscore, double countDP, double fdr) {
        this.strpeptide = strpeptide;
        this.msonefeature = msonefeature;
        this.proteinname = proteinname;
        this.decoyTarget = decoyTarget;
        this.pepmass = pepmass;
        this.msonemass = msonemass;
        this.msonetime = msonetime;
        this.msoneMZ = msoneMZ;
        this.msoneZ = msoneZ;
        this.BIon = BIon;
        this.YIon = YIon;
        this.targetscore = targetscore;
        this.countDP = countDP;
        this.fdr = fdr;
    }

    @Override
    public String toString() {
        return
                 strpeptide +
                "\t" + msonefeature +
                "\t" + proteinname +
                "\t" + decoyTarget +
                "\t" + pepmass +
                "\t" + msonemass +
                "\t" + msonetime +
                "\t" + msoneMZ +
                "\t" + msoneZ +
                "\t" + BIon +
                "\t" + YIon +
                "\t" + targetscore +
                "\t" + countDP +
                "\t" +fdr +
                "\t" + sharepeaksMax+
                 "\t" + sharePeaksAll+
                "\t" + strSharePeaksMaxInfo+
                "\t" + strSharePeaksAllInfo
                ;

    }



    @Override
    public int hashCode() {
        return Objects.hash(msonefeature);
    }

    @Override
    public int compareTo(PSMResult o) {
        return this.msonefeature.compareTo(o.msonefeature);
    }
}
