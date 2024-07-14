package edu.uw.waterlooms.match;

import java.util.ArrayList;


public class Genome {
    public String id;
    public String description;
    public ArrayList<Peptide> composition;

    // used mainly to do tryptic digestion
    public Genome(String id, String composition, String description) {
        this.id = id;
        this.description = description;
        this.composition = TrypsinRule(id, composition);
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
            if (cur == 'K' || cur == 'R') {
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
}
