package edu.uw.waterlooms.peptideMatch;



import java.util.ArrayList;


// stores the score of a sliding window and peaks matched in it
public class MatchedInterval implements Comparable<MatchedInterval>{
    Double score;
    ArrayList<Peak> peaks;

    public MatchedInterval(double score, ArrayList<Peak> peaks) {
        this.score = score;
        this.peaks = peaks;
    }

    @Override
    public int compareTo(MatchedInterval other) {
        Double score = this.score;
        return score.compareTo(other.score);
    }

}
