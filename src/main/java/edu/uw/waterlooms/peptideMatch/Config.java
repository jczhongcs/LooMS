package edu.uw.waterlooms.peptideMatch;


public class Config {
    static Enums.RunMode mode = Enums.RunMode.DEBUG;
    public static Enums.ScoreMode scoreMode = Enums.ScoreMode.RELATIVE;
    public static String userHome = System.getProperty("user.home");
//    public static String fastaFolder = "/tmp/data/mouse/";
    public static String fastaFolder = "/tmp/data/";
//    public static String spectrumFolder = "\\Desktop\\UWAT\\fourth_year\\summer\\spec_filter\\data\\spectrums\\";
//    public static String spectrumFolder = "/tmp/data/mouse/";//MOUSE DATASET
    public static String spectrumFolder = "/tmp/data/";//HUMAN DATASET
    public static String spectrumFolderLocal = "/Users/jianzhong/Documents/uwaterloo_dia/DIA Data/";

    //    public static double maxGap = 0;
    public static double tsThreshold = 0.33333333;  // minutes
}
