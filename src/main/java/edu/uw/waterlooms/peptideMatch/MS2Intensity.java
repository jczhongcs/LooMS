package edu.uw.waterlooms.peptideMatch;


import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class MS2Intensity {

    String strPep;
    int iLen;
    double dRT;
    public double dPredictRT;
    int icharge;
    public double[][]  arrdIntensity;

    public int[] arrbTop6IntensityPosCharge1;
    public int[] arryTop6IntensityPosCharge1;
    //    public int[] arrallTop6IntensityCharge1;
    public int arrallTop6IntensityCharge1;
    public int[] arrallTop6IntensityCharge1Pos;

    public int[] arrallTop6IntensityCount;
    public int[] arrallTop6IntensityPos;
    public void generatePosforTopN()
    {
        int k = 6;

        double[] arrb1ion = new double[arrdIntensity.length];
        double[] arry1ion = new double[arrdIntensity.length];
        for(int i =0;i<arrdIntensity.length;i++)
        {
            arrb1ion[i] = arrdIntensity[i][0];
            arry1ion[i] = arrdIntensity[i][2];
        }

        if(arrdIntensity.length>=k)
        {
            arrbTop6IntensityPosCharge1 = indexesOfTopElements(arrb1ion,k);//b
            arryTop6IntensityPosCharge1 = indexesOfTopElements(arry1ion,k);//y

        }else
        {
            arrbTop6IntensityPosCharge1 = new int [arrdIntensity.length];

            for(int i =0;i<arrdIntensity.length;i++)
            {
                arrbTop6IntensityPosCharge1[i] = i;
            }
            arryTop6IntensityPosCharge1 = new int [arrdIntensity.length];

            for(int i =0;i<arrdIntensity.length;i++)
            {
                arryTop6IntensityPosCharge1[i] = i;
            }
        }
        double arrCombine[] = new double[arrb1ion.length+arry1ion.length];
        for(int i=0;i<arrb1ion.length;i++)
        {
            arrCombine[i]=arrb1ion[i];
        }
        for(int i=0;i<arry1ion.length;i++)
        {
            arrCombine[i+arrb1ion.length]=arry1ion[i];
        }
        int[] indicesCombine = indexesOfTopElements(arrCombine,k);
//        arrallTop6IntensityCharge1 = new int[k];
        arrallTop6IntensityCharge1 = 0;//new int[k];
        arrallTop6IntensityCharge1Pos = new int[k];
        for(int i = 0; i < indicesCombine.length; i++) {
            if(indicesCombine[i]<arrb1ion.length)
            {
                arrallTop6IntensityCharge1++;
//                arrallTop6IntensityCharge1[i]=0;
//                arrallTop6IntensityPos[i] = arrbTop6IntensityPosCharge1[indicesCombine[i]];
                arrallTop6IntensityCharge1Pos[i] = indicesCombine[i];
            }else
            {
//                arrallTop6IntensityCharge1[i]=2;
                arrallTop6IntensityCharge1Pos[i] =  indicesCombine[i]  - arrb1ion.length;
            }
        }


        double arrallCombine[] = new double[4*arrb1ion.length];
        for(int i=0;i<arrdIntensity.length;i++)
        {
            for(int j =0;j<4;j++)
            {
                arrallCombine[i+j*arrdIntensity.length]=arrdIntensity[i][j];
            }
        }

        int[] indicesallCombine = indexesOfTopElements(arrallCombine,k);
//        arrallTop6IntensityCharge1 = new int[k];
        arrallTop6IntensityCount = new int[4];//new int[k];
        arrallTop6IntensityPos = new int[k];
        int j=1;
        for(int i = 0; i < indicesallCombine.length; i++) {
            if(indicesallCombine[i]<arrb1ion.length*j)
            {
                arrallTop6IntensityCount[j-1]++;
//                arrallTop6IntensityCharge1[i]=0;
//                arrallTop6IntensityPos[i] = arrbTop6IntensityPosCharge1[indicesCombine[i]];
                arrallTop6IntensityPos[i] = indicesallCombine[i] - arrb1ion.length*(j-1);
            }else
            {
                arrallTop6IntensityCount[j]++;

//                arrallTop6IntensityCharge1[i]=2;
                arrallTop6IntensityPos[i] =  indicesallCombine[i]  - arrb1ion.length*j;
                j++;
            }
        }


    }
    public static void main(String[] args) {
        double[][] arr = {
                {0.9,      3.2,        23.12,      7.2},
                {0.9,      112.0,      23.12,      0.86},
                {3.09,      39,         223.12,      120},
                {0.9,      3.2,        23.12,      7.2},
                {0.9,      112.0,      23.12,      0.86}
                ,
                {1.227,        39,         68.81,       10},
                {3.09,      3.2,        5.2,        7.2},
                {3.09,      203.12,     10.4,       0.86},
                {3.09,      9,          18,         10}
        };
        MS2Intensity mi = new MS2Intensity();
        mi.arrdIntensity = arr;
        mi.generatePosforTopN();
        System.out.println("te");
//        int[] indexes = indexesOfTopElements(arr,6);
        int[] indexes = mi.arrallTop6IntensityCharge1Pos;
        for(int i = 0; i < indexes.length; i++) {
            int index = indexes[i];
            System.out.println(index + " " + mi.arrdIntensity[index][i<mi.arrallTop6IntensityCharge1?0:2]);
        }
    }
    //READ pDeep3 data with specical format like multiple 空格
    public static Map<String,MS2Intensity> getMapMS2IntensityFromProsit(String filepath) throws IOException {
//    public static void main(String[] args) throws IOException {
        // write your code here
        System.out.println("prositstart");

        long start = System.currentTimeMillis();


//        Map<Integer,MS2Intensity> mapStrMS2Intensity = new HashMap<>();
        Map<String,MS2Intensity> mapStrMS2Intensity = new HashMap<>();
        FileReader freader;
        BufferedReader br;
//        FileReader freaderConcise;
//        BufferedReader brConcise;
        long timeSpan = 0;

//        System.out.println(filepath);
        List<File> filesInFolder = Files.walk(Paths.get(filepath))
                .filter(Files::isRegularFile)
                .map(Path::toFile)
                .collect(Collectors.toList());

//        for (String strPrositResultFileName:Utils.strPrositResultFileName) {
        for (File strPrositResultFileName:filesInFolder) {
//            String dataFile = Config.spectrumFolder + strPrositResultFileName; // Utils.strpDeepOutfile;//pdeep3output.txt";
//            String dataFile = filepath+strPrositResultFileName; // Utils.strpDeepOutfile;//pdeep3output.txt";
//        String dataFile = Config.spectrumFolder +Utils.strPrositResultFileName; // Utils.strpDeepOutfile;//pdeep3output.txt";
//        String dataFile = Config.spectrumFolder + Utils.strpDeepOutfile;//+"result-filename"; // Utils.strpDeepOutfile;//pdeep3output.txt";
//        String conciseFile = Config.spectrumFolder +"pDeep3/newVersion/pDeep3/pDeep/tmp/predict/peptide.txt"; //"pDeep3/newVersion/pDeep3/peptidesForpProsit_0502.csv";//"alloutputconcise.csv";

//        String strTest="AMGDQMLQVTSSTLGDVLETENSNYARKEEEAELR3";
            String line;
            int i = 1;

            try {
                freader = new FileReader(strPrositResultFileName);
                br = new BufferedReader(freader);
//            freaderConcise = new FileReader(conciseFile);
//            brConcise = new BufferedReader(freaderConcise);
//            brConcise.readLine();
                String strLine;


                while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) ;//定位到当前位置

                while (line != null) {
//                String[] strPepInfo = strLine.split("\t");


                    String[] strPredictPepInfo = line.split("\\|");
                    double[][] arrPredictMSMS = new double[strPredictPepInfo[1].length() - 1][4];//length -1 with b b+2 y y+2
                    StringBuffer buffer = new StringBuffer();
//                if(!strPepInfo[0].equals(strPredictPepInfo[1]))
//                {
//                    System.out.println(strPepInfo[0]+"    "+strPredictPepInfo[1]);
//                }
                    while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) {

                        int nextPos = line.indexOf(" ", 0);
                        String strDoubleIntensity = line.substring(0, nextPos);
                        String iontype = line.substring(nextPos + 1, nextPos + 2);

                        int nextPosPlus = line.indexOf(" ", nextPos + 1);
                        int icharge = 1;
                        if (nextPosPlus > 0) {
                            String strCharge = line.substring(nextPosPlus + 1, nextPosPlus + 2);
                            if (strCharge.length() > 0)
                                icharge = Integer.parseInt(strCharge);
                            else
                                System.out.println("debug");
                        } else {
                            nextPosPlus = line.length();
                        }
                        String strIonPos = line.substring(nextPos + 2, nextPosPlus);

                        int iIonPos = Integer.parseInt(strIonPos);
                        if (icharge <= 2) {
                            if (iontype.equals("b")) {
                                icharge--;
                            } else {
                                icharge++;

                            }
                            //过滤3charge data
                            arrPredictMSMS[iIonPos - 1][icharge] = Double.parseDouble(strDoubleIntensity);

                        }

                    }


                    MS2Intensity ms2Inten = new MS2Intensity();
                    ms2Inten.dPredictRT = Double.parseDouble(strPredictPepInfo[3].substring(1, strPredictPepInfo[3].length() - 2));
                    ms2Inten.arrdIntensity = arrPredictMSMS;
                    ms2Inten.generatePosforTopN();
//                mapStrMS2Intensity.put(strPepInfo[0]+strPepInfo[2], ms2Inten);
                    mapStrMS2Intensity.put(strPredictPepInfo[1] + strPredictPepInfo[2], ms2Inten);
                    i++;
                    if (i % 1000000 == 0) {
                        timeSpan = System.currentTimeMillis() - start;
                        System.out.println("Time:" + timeSpan + " i=" + i);
                    }
//                if(line==null) break;
                }

            } catch (FileNotFoundException noFile) {
                throw new FileNotFoundException();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }


            br.close();
            freader.close();
        }
        timeSpan = System.currentTimeMillis() - start;
        System.out.println(timeSpan);
        return mapStrMS2Intensity;
//        MS2Intensity test = mapStrMS2Intensity.get(strTest);
//        System.exit(0);

    }
    public static Map<String,MS2Intensity> getMapMS2IntensityFromProsit() throws IOException {
//    public static void main(String[] args) throws IOException {
        // write your code here
        System.out.println("prositstart");

        long start = System.currentTimeMillis();


//        Map<Integer,MS2Intensity> mapStrMS2Intensity = new HashMap<>();
        Map<String,MS2Intensity> mapStrMS2Intensity = new HashMap<>();
        FileReader freader;
        BufferedReader br;
//        FileReader freaderConcise;
//        BufferedReader brConcise;
        long timeSpan = 0;

        for (String strPrositResultFileName:Utils.strPrositResultFileName) {
//            String dataFile = Config.spectrumFolder + strPrositResultFileName; // Utils.strpDeepOutfile;//pdeep3output.txt";
            String dataFile = strPrositResultFileName; // Utils.strpDeepOutfile;//pdeep3output.txt";
//        String dataFile = Config.spectrumFolder +Utils.strPrositResultFileName; // Utils.strpDeepOutfile;//pdeep3output.txt";
//        String dataFile = Config.spectrumFolder + Utils.strpDeepOutfile;//+"result-filename"; // Utils.strpDeepOutfile;//pdeep3output.txt";
//        String conciseFile = Config.spectrumFolder +"pDeep3/newVersion/pDeep3/pDeep/tmp/predict/peptide.txt"; //"pDeep3/newVersion/pDeep3/peptidesForpProsit_0502.csv";//"alloutputconcise.csv";

//        String strTest="AMGDQMLQVTSSTLGDVLETENSNYARKEEEAELR3";
            String line;
            int i = 1;

            try {
                freader = new FileReader(dataFile);
                br = new BufferedReader(freader);
//            freaderConcise = new FileReader(conciseFile);
//            brConcise = new BufferedReader(freaderConcise);
//            brConcise.readLine();
                String strLine;


                while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) ;//定位到当前位置

                while (line != null) {
//                String[] strPepInfo = strLine.split("\t");


                    String[] strPredictPepInfo = line.split("\\|");
                    double[][] arrPredictMSMS = new double[strPredictPepInfo[1].length() - 1][4];//length -1 with b b+2 y y+2
                    StringBuffer buffer = new StringBuffer();
//                if(!strPepInfo[0].equals(strPredictPepInfo[1]))
//                {
//                    System.out.println(strPepInfo[0]+"    "+strPredictPepInfo[1]);
//                }
                    while ((line = br.readLine()) != null && !(line.startsWith(">peptide|"))) {

                        int nextPos = line.indexOf(" ", 0);
                        String strDoubleIntensity = line.substring(0, nextPos);
                        String iontype = line.substring(nextPos + 1, nextPos + 2);

                        int nextPosPlus = line.indexOf(" ", nextPos + 1);
                        int icharge = 1;
                        if (nextPosPlus > 0) {
                            String strCharge = line.substring(nextPosPlus + 1, nextPosPlus + 2);
                            if (strCharge.length() > 0)
                                icharge = Integer.parseInt(strCharge);
                            else
                                System.out.println("debug");
                        } else {
                            nextPosPlus = line.length();
                        }
                        String strIonPos = line.substring(nextPos + 2, nextPosPlus);

                        int iIonPos = Integer.parseInt(strIonPos);
                        if (icharge <= 2) {
                            if (iontype.equals("b")) {
                                icharge--;
                            } else {
                                icharge++;

                            }
                            //过滤3charge data
                            arrPredictMSMS[iIonPos - 1][icharge] = Double.parseDouble(strDoubleIntensity);

                        }

                    }


                    MS2Intensity ms2Inten = new MS2Intensity();
                    ms2Inten.dPredictRT = Double.parseDouble(strPredictPepInfo[3].substring(1, strPredictPepInfo[3].length() - 2));
                    ms2Inten.arrdIntensity = arrPredictMSMS;
                    ms2Inten.generatePosforTopN();
//                mapStrMS2Intensity.put(strPepInfo[0]+strPepInfo[2], ms2Inten);
                    mapStrMS2Intensity.put(strPredictPepInfo[1] + strPredictPepInfo[2], ms2Inten);
                    i++;
                    if (i % 1000000 == 0) {
                        timeSpan = System.currentTimeMillis() - start;
                        System.out.println("Time:" + timeSpan + " i=" + i);
                    }
//                if(line==null) break;
                }

            } catch (FileNotFoundException noFile) {
                throw new FileNotFoundException();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }


            br.close();
            freader.close();
        }
        timeSpan = System.currentTimeMillis() - start;
        System.out.println(timeSpan);
        return mapStrMS2Intensity;
//        MS2Intensity test = mapStrMS2Intensity.get(strTest);
//        System.exit(0);

    }
    static int[] indexesOfTopElements(double[] orig, int nummax) {
        double[] copy = Arrays.copyOf(orig,orig.length);
        Arrays.sort(copy);

        double dminValue=copy[copy.length- nummax];
        int iCountleastSameCount = 1;
        for(int i=copy.length- nummax+1;i<copy.length && dminValue==copy[i];i++,iCountleastSameCount++);

        int[] result = new int[nummax];
        int resultPos = 0;
        for(int i = 0; i < orig.length; i++) {
            double onTrial = orig[i];
            if(onTrial > dminValue ) {

                result[resultPos++] = i;

            }else  if(onTrial == dminValue){
                if(iCountleastSameCount > 0)
                {
                    result[resultPos++] = i;
                    iCountleastSameCount--;
                }
            }

        }

        return result;
    }


}
