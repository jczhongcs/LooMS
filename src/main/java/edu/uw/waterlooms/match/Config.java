package edu.uw.waterlooms.match;

public class Config {
    static Enums.RunMode mode = Enums.RunMode.DEBUG;
    static Enums.ScoreMode scoreMode = Enums.ScoreMode.RELATIVE;
    public static String userHome = System.getProperty("user.home");
    public static String spectrumFolder = "\\Desktop\\UWAT\\fourth_year\\summer\\spec_filter\\data\\spectrums\\sample\\";

    public static String testspectrumFolder = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/";


    //public static double tsThreshold = 0.33333333; // minutes
    // Sliding RT timestamp value in MINUTES ~ 9.6 seconds
    public static double tsThreshold = 0.1666666;
    //public static double tsThreshold = 0.08333333;

}
