package edu.uw.waterlooms.peptideMatch;


import java.util.ArrayList;
import java.util.List;
import java.util.TreeMap;
import java.util.stream.Collectors;


public class Genome {
    public String id;
    public String description;
    public ArrayList<Peptide> arrPeps;

    // used mainly to do tryptic digestion
    public Genome(String id, String arrPeps, String description) {
        this.id = id;
        this.description = description;
        if(Utils.bTrypsinRule) {
            this.arrPeps = TrypsinRule(id, arrPeps);
        }
        else if(Utils.bTrypsinRuleWithoutP) {
            this.arrPeps = TrypsinRuleNoConsiderP(id, arrPeps);
        }
        if(Utils.bzGenerateSemiTrypsin)
        {
            semiTrypsinWithNTerminate();
            semiTrypsinWithCTerminate();
        }
        if(Utils.bzGenerateMissCleavage)
        {
            missCleavageTrypsin();
        }

    }
    public static ArrayList<Peptide> TrypsinRuleNoConsiderP(String id, String composition) {
        ArrayList<Peptide> res = new ArrayList<>();
        int len = composition.length();

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < len; ++i) {
            char cur = composition.charAt(i);
            // check if is end
            if (cur == '*') {
                // do not append empty peptides
                if(sb.length() > 0) {
                    res.add(new Peptide(id, sb.toString()));
                    sb.setLength(0);
                }
                break;
            }

            // check illegal params; if seen illegal character, cut it out to make 2 peptides
            if (!AminoAcid.mass.containsKey(cur)) {
                // do not append empty peptides
                if(sb.length() == 0) {
                    continue;
                }
                res.add(new Peptide(id, sb.toString()));
                sb.setLength(0);
                continue;
            }

            sb.append(cur);
            if (cur == 'K' || cur == 'R') {//[KR]
                // trypsin digestion
//                if (i < len -1) {
//                    if (composition.charAt(i+1) != 'P') {
                        res.add(new Peptide(id, sb.toString()));
                        sb.setLength(0);
//                    }
//                }
            }
        }

        // last one
        if(sb.length() > 0) {
            res.add(new Peptide(id, sb.toString()));
            sb.setLength(0);
        }

        return res;
    }
    public static ArrayList<Peptide> TrypsinRule(String id, String composition) {
        ArrayList<Peptide> res = new ArrayList<>();
        int len = composition.length();

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < len; ++i) {
            char cur = composition.charAt(i);
            // check if is end
            if (cur == '*') {
                // do not append empty peptides
                if(sb.length() > 0) {
                    res.add(new Peptide(id, sb.toString()));
                    sb.setLength(0);
                }
                break;
            }

            // check illegal params; if seen illegal character, cut it out to make 2 peptides
            if (!AminoAcid.mass.containsKey(cur)) {
                // do not append empty peptides
                if(sb.length() == 0) {
                    continue;
                }
                res.add(new Peptide(id, sb.toString()));
                sb.setLength(0);
                continue;
            }

            sb.append(cur);
            if (cur == 'K' || cur == 'R') {//[KR] NOT [RP or KP]
                // trypsin digestion
                if (i < len -1) {
                    if (composition.charAt(i+1) != 'P') {
                        res.add(new Peptide(id, sb.toString()));
                        sb.setLength(0);
                    }
                }
            }
        }

        // last one
        if(sb.length() > 0) {
            res.add(new Peptide(id, sb.toString()));
            sb.setLength(0);
        }

        return res;
    }
    public void semiTrypsinWithNTerminate()
    {
        ArrayList<Peptide> semiTrypsinWithNTerminatePeptides = new ArrayList<>();
        for(Peptide pep: arrPeps)
        {
            for (int ilen = 1; (pep.composition.length()-ilen)>=Utils.thresholdPeptideSizeMin; ilen++)
            {
                semiTrypsinWithNTerminatePeptides.add(new Peptide(pep.id, pep.composition.substring(ilen),1.0-((double) (ilen)/pep.composition.length())));
            }
        }
        arrPeps.addAll(semiTrypsinWithNTerminatePeptides);


    }
    public void semiTrypsinWithCTerminate()
    {
        ArrayList<Peptide> semiTrypsinWithCTerminatePeptides = new ArrayList<>();
        for(Peptide pep: arrPeps)
        {
            for (int ilen = 1; (pep.composition.length()-ilen)>=Utils.thresholdPeptideSizeMin; ilen++)
            {
                semiTrypsinWithCTerminatePeptides.add(new Peptide(pep.id, pep.composition.substring(0,pep.composition.length()-ilen),1.0-((double) (ilen)/pep.composition.length())));
            }
        }
        arrPeps.addAll(semiTrypsinWithCTerminatePeptides);
    }
    public void missCleavageTrypsin()
    {
        ArrayList<Peptide> missCleavageTrypsinPeptides = new ArrayList<>();
        for(int i = 0;i<arrPeps.size();i++)
        {
            String pepComposition = arrPeps.get(i).composition;
            int iJoin = 0;//miss cut count

            for(int iCleav = 0;iCleav<Utils.iCleavageThreshold;iCleav++)
            {
                if(( i+iCleav+1<arrPeps.size()) && arrPeps.get(i).id.equals(arrPeps.get(i+iCleav+1).id))//是否是同一个protein
                {
                    if((pepComposition.length() + arrPeps.get(i+iCleav+1).composition.length()) < Utils.thresholdPeptideSizeMax)
                    {
                        pepComposition += arrPeps.get(i+iCleav+1).composition;
                        iJoin++;
                    }else
                    {
                        break;
                    }

                }

                if( iJoin > 0 )
                {
                    missCleavageTrypsinPeptides.add(new Peptide(arrPeps.get(i).id, pepComposition,1.0-((double) (iJoin)/pepComposition.length())));

                }
            }



        }

        arrPeps.addAll(missCleavageTrypsinPeptides);
    }
//    public void removePeptideFromBothTargetAndDecoy()
//    {
//        TreeMap<String, List<Peptide>> mapPepstrListPep = arrPeps.stream().sorted()
//                .collect(Collectors.groupingBy(Peptide::getComposition, TreeMap::new, Collectors.toList()));
//
//        ArrayList<Peptide> arrPepsRemove = new ArrayList<>();
//        for ( String strPep: mapPepstrListPep.keySet()) {
//
//
//            boolean bsameType = false;
//            if(mapPepstrListPep.get(strPep).size()>1) {
//                int i = 0;
//                //用异或来判断peptide是否来自于target和decoy
//                for (Peptide pep : mapPepstrListPep.get(strPep)) {
//
//                    if (i == 0) {
//                        bsameType = pep.composition.contains("DeBruijn");
//                        i++;
//                    } else {
//                        bsameType ^= pep.composition.contains("DeBruijn");
//                    }
//                }
//            }
//            if(!bsameType) arrPepsRemove.addAll(mapPepstrListPep.get(strPep));
//
//        }
//        arrPeps = arrPepsRemove;
//    }

}
