package edu.uw.waterlooms.entity;

import java.util.List;
import java.util.Map;

public class SignalIsotope {

    public int iRTpos;
    public int iMzPos;
    public double mzvalue;

    public List<Integer> liNumber_Isotope;
    public Map<Integer,List<Integer>> lliIstopePos;
    public Map<Integer,List<Double>> lliIstopeValue;
}
