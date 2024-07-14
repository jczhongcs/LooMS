package edu.uw.waterlooms.match;

import edu.uw.waterlooms.entity.XIC;

import java.util.ArrayList;
import java.util.List;


// stores the score of a sliding window and peaks matched in it
public class MatchedInterval implements Comparable<MatchedInterval>{
    MatchFeatures matchedFeatures;

    Double score;
    List<XIC> peaks;
    XIC precursor;


    public MatchedInterval(MatchFeatures matchedFeatures){
        this.matchedFeatures = matchedFeatures;
    }

    public MatchedInterval(double score, List<XIC> peaks) {
        this.score = score;
        this.peaks = peaks;
        this.precursor = null;
    }

    public MatchedInterval(double score, List<XIC> peaks, XIC precursor){
        this.score = score;
        this.peaks = peaks;
        this.precursor = precursor;
    }

    public XIC returnFirstTrail(){
        // TODO: Fix this code smell
        return this.peaks.get(0);
    }

    public MatchFeatures getMatchedFeatures(){
        return this.matchedFeatures;
    }

    public double getScore(){
        return this.score;
    }

    public void setScore(Double score) {
        this.score = score;
    }

    public XIC getPrecursor() { return precursor; }

    public void setPrecursor(XIC precursor) { this.precursor = precursor; }

    @Override
    public int compareTo(MatchedInterval other) {
        Double score = this.score;
        return score.compareTo(other.score);
    }

}
