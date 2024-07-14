package edu.uw.waterlooms.entity;

import edu.uw.waterlooms.match.Utils;
import edu.uw.waterlooms.neuralscoring.RetentionTimeMetadataTuple;
import org.apache.commons.math3.util.Pair;

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;

public class IsolationWindow implements Serializable {
    // map of retention time (rt) to max intensity at this rt within this window
    public HashMap<Double, Double> max; // = new HashMap<>();
    // map of RT to Intensities Observed at the RT (used to compute RANK)
    public HashMap<Double, ArrayList<Double>> rtToIntensitySet;

    // map of unique RTs, along with the XIC/TRAILs that fall within that RT +/- some offset (default 20s) for matching
    public HashMap<Double, ArrayList<XIC>> rtToTrailSet;

    // map of unique RTs, containing all the precursor(s) that fall within RT - some offset (default 10s) for precursros
    public HashMap<Double, ArrayList<XIC>> rtToPrecursorSet;

    /* FOR NEURAL SCORING */
    public HashMap<XIC, Integer> trailToMatchedPeptidesSet;

    // map of unique RTs, to some metadata regarding the peaks falling within that RT; used for Neural Scoring
    public HashMap<Double, RetentionTimeMetadataTuple> rtToMetadataSet;

    // trails are to be sorted by mass
    public ArrayList<XIC> xics; // = new ArrayList<>();
    public double mzLow;
    public double mzHigh;
    public int iNo;

    public IsolationWindow(double mzLow, double mzHigh)
    {
        this.xics = new ArrayList<>();
        this.max = new HashMap<>();
        this.rtToIntensitySet = new HashMap<>();
        this.rtToTrailSet = new HashMap<>();
        this.rtToPrecursorSet = new HashMap<>();
        this.trailToMatchedPeptidesSet = new HashMap<>();
        this.rtToMetadataSet = new HashMap<>();
        this.mzLow = mzLow;
        this.mzHigh = mzHigh;
    }

    public IsolationWindow(double mzLow, double mzHigh,int NO)
    {
        this(mzLow,mzHigh);
        this.iNo = NO;
    }
    // given mass, returns an array of indices of all matched trails
    public ArrayList<Integer> FindTrails(double mass, int charge) {
        ArrayList<Integer> res = new ArrayList<>();
        double mz = Utils.MassToMz(mass, charge);
        XIC xic = new XIC(mz);

        int index = Collections.binarySearch(xics, xic);

        if (index >= 0) {
            int first = index;
            int last = index;

            res.add(index);

            while (first > 0 && xics.get(first - 1).compareTo(xic) == 0) {
                res.add(first - 1);
                --first;
            }
            while (last < xics.size()-1 && xics.get(last + 1).compareTo(xic) == 0) {
                res.add(last + 1);
                ++last;
            }

            return res;
        }
        return null;
    }

    // given rt's and intensities of a new trail, update max
    public void PutRtMax(Double[] rts, Double[] ints) {
        for (int i = 0; i < rts.length; ++i) {
//            Double key = Utils.RoundToX(rts[i], Config.maxGap);
            Double key = rts[i];
            Double val = ints[i];
            if (max.containsKey(key)) {
                if (max.get(key) < val) {
                    max.put(key, val);
                }
            } else {
                max.put(key, val);
            }
        }
    }

    /**
     * Return a map of the retention time to the max intensity of trial at that RT
     * @return
     */
    public Map<Double, Double> rtToMaxIntensityMap(){
        return xics.stream().collect(Collectors.toMap(
           XIC::getRtAtMaxIntensity, XIC::getMaxIntensity
        ));
    };

    /**
     * Utility method to populate the RT : Precursor hashmap.
     * For each RT, sets the possible precursors in the - toleranceInMinutes given the RT
     * @param toleranceInMinutes
     * @param precursors
     */
    public void populatertToPrecursorSet(double toleranceInMinutes, ArrayList<XIC> precursors){

        // Go through the list of filtered Trails
        this.rtToTrailSet.forEach((rt, xic_arraylist) -> {
            // Add these RTs to the precursor set for precursor subsetting
            if(!this.rtToPrecursorSet.containsKey(rt)){
                XIC upperRTXIC = new XIC(0d, rt);
                XIC lowerRTXIC = new XIC(0d, rt - toleranceInMinutes);
                int upperIndex = parseInsertionPoint(Collections.binarySearch(precursors, upperRTXIC, XIC::compareRtXIC));
                int lowerIndex = parseInsertionPoint(Collections.binarySearch(precursors, lowerRTXIC, XIC::compareRtXIC));
                this.rtToPrecursorSet.put(rt, new ArrayList<XIC>(precursors.subList(lowerIndex, upperIndex)));
            }
        });

        // Once all Precursors have been selected by RT - tolerance, sort the resulting precursor set by M/Z
        this.rtToPrecursorSet.forEach((rt, trails) -> {
            // Sort
            trails.sort(XIC::compareMzXIC);
        });
    }

    /**
     * Utility method to populate the RT : [TRAIL] hashmap.
     * @param toleranceInMinutes int of seconds of the tolerance around the RT to be subsetted
     */
    public void populatertToTrailSet(double toleranceInMinutes){
        // Sort the internal XICs by RT
        this.xics.sort(XIC::compareRtXIC);

        /*
         * This may be quite MEMORY inefficient
         * For each XIC/Trail ...
         * Read the RT at max intensity
         * IF this RT is at least 1/2 the tolerance away from the previously seen RT...
         * If this is a unique RT, add the trail to the hashmap
         * After the additions, perform a sort on each one of the element(s)
         * Use this rtToTrailSet for scoring purposes as you no longer need to sort
         */
        AtomicReference<Double> seenRT = new AtomicReference<>((double) 0);
        this.xics.forEach(xic ->{
            double currentRT = xic.getRtAtMaxIntensity();

            // Skip if the RT differential (half the given tolerance) isn't great enough
            if (Math.abs(currentRT - seenRT.get()) < (toleranceInMinutes/2d)){
                return;
            }
            seenRT.set(currentRT);
            // If the RT has not been observed before ...
            if(!this.rtToTrailSet.containsKey(currentRT)){

                // Subset all RTs around the secondTolerance
                ArrayList<XIC> subsettedTrails = subsetXICBasedOnRT(xic, toleranceInMinutes);
                this.rtToTrailSet.put(currentRT, subsettedTrails);
            }
            // Otherwise skip it, we only want the unique RTs
        });

        // Once all subsets have been sufficiently computed, sort the resultant XIC subsets.
        this.rtToTrailSet.forEach((rt, trails) -> {
            // Sort
            trails.sort(XIC::compareMzXIC);
            // TODO: See if you need to re-set the key'd values
        });
    }

    // TODO: May need to debug this issue
    public void populateMax(){
        // For each XIC belonging to the collection of XICs
        this.xics.forEach(xic ->
        {
            // For EACH retention time within the XIC
            ArrayList rts = xic.getRetentionTimes();
            ArrayList ints = xic.getIntensities();
            for (int rtIndex = 0; rtIndex < xic.getRetentionTimes().size(); rtIndex ++)
            {
                Double tmpRT = (double) rts.get(rtIndex);
                Double tmpIT = (double) ints.get(rtIndex);

                // If the RT does not exist
                if (!this.max.containsKey(tmpRT))
                {
                    // Add the RT and Intensity (RetentionTime, Intensity)
                    this.max.put(tmpRT, tmpIT);
                    this.rtToIntensitySet.put(tmpRT,
                            new ArrayList<>()); // 1st such occurrence for tmpRT
                } else {
                    // IF the RT Exists
                    // Append the new Intensity to the set
                    this.rtToIntensitySet.get(tmpRT).add(tmpIT);

                    // Update running sum for Intensity at the given RT
                    //double currIntensitySum = this.max.get(tmpRT);
                    //this.max.replace(tmpRT, tmpIT); // Replace with a higher intensity

                    // IF the RT exists, and the current Intensity is LOWER
                    if (this.max.get(tmpRT) < tmpIT)
                    {
                        // Replace with max Intensity
                        this.max.replace(tmpRT, tmpIT);
                    }
                }
            }

            /*
            for (int rtIndex = 0; rtIndex < xic.getRetentionTimes().size(); rtIndex ++){
                Double tmpRT = (double) rts.get(rtIndex);
                Double tmpIT = (double) ints.get(rtIndex);
                // Normalize each value
                double runningSum = this.max.get(tmpRT);
                int numSeen = this.rtCount.get(tmpRT);

                // Normalize the SUM of the Intensities at a given RT, with the number of peaks belonging to that RT
                this.max.replace(tmpRT, (runningSum / (double) numSeen));
            }*/

        });

        // For each KEY in the hashmap, sort the value(s) by intensity
        // This is for
        Iterator it = this.rtToIntensitySet.entrySet().iterator();
        while (it.hasNext()){
            Map.Entry pair = (Map.Entry)it.next();
            double rt = (double)pair.getKey();
            ArrayList<Double> intensitiesAtGivenRT = this.rtToIntensitySet.get(rt);
            Collections.sort(intensitiesAtGivenRT);
            // Reset the intensity set
            this.rtToIntensitySet.replace(rt, intensitiesAtGivenRT);
        }

        /*
        // TODO: This normalization does not work.
        // TODO: This will result in a SMALLER maxIntensity than a peakIntensity
        // For each KEY in the MAX hashmap (keys are the Retention Times observed)
        Iterator it = this.max.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();

            // Get the Retention Time
            double rt = (double) pair.getKey();
            // Get the number of peaks contributing to this intensity summation
            int numSeen = this.rtCount.get(rt);
            // Get the running intensity sunmmation
            double runningsum = this.max.get(rt);

            // Normalize it to the average intensity at the given RT
            this.max.replace(rt, (runningsum/(double)numSeen));
            //it.remove(); // avoids a ConcurrentModificationException
        }
        */
    }

    public void removeDuplicateTrails(){
        // For each xic entry in this object
        // Ensure that there are no duplicate(s)

        int y = 0;

        // First index is the RT
        // Second index is the MZ
        // If the mz is longer, then remove the set
        // At the end collect all to a single set
        HashMap<Double, HashMap<Double, XIC>> rt_mz_xic_map = new HashMap<>();

        this.xics.forEach(
            xic -> {
                // RT does not exist, add to hashmap & create new hashmap
                if(!rt_mz_xic_map.containsKey(xic.getRtAtMaxIntensity()))
                {
                    // Create the RT KEY, and MAP double
                    HashMap<Double, XIC> mzToXICMap = new HashMap<>();
                    mzToXICMap.put(xic.getMZAtMaxIntensity(), xic);
                    rt_mz_xic_map.put(
                            xic.getRtAtMaxIntensity(),
                            mzToXICMap);
                } else {
                    // RT Exists, Check MZ
                    HashMap<Double, XIC> secondMap = rt_mz_xic_map.get(xic.getRtAtMaxIntensity());
                    // New [RT [MZ, XIC]] entry
                    if (!secondMap.containsKey(xic.getMZAtMaxIntensity())){
                        // MZ Doesn't exist, create new XIC
                        secondMap.put(xic.getMZAtMaxIntensity(), xic);
                    } else {
                        // RT exists, and MZ exists, compare the length of the XICs
                        if (
                                // Length of the current xic is smaller than the proposed one
                                secondMap.get(xic.getMZAtMaxIntensity()).getRetentionTimes().size() <
                                        xic.getRetentionTimes().size()
                        ){
                            secondMap.replace(xic.getMZAtMaxIntensity(), xic);
                        }
                    }


                }
            }
        );


        ArrayList<XIC> deduplicatedXICs = new ArrayList<>();

        // Iterate through and collect the "deduplicated" set
        for(Map.Entry<Double, HashMap<Double, XIC>> rtkey : rt_mz_xic_map.entrySet()){
            HashMap<Double,XIC> mzkey = rtkey.getValue();
            for(Map.Entry<Double, XIC> filteredXICs : mzkey.entrySet()){
                deduplicatedXICs.add(filteredXICs.getValue());
            }
        }

        // Assign the deduplicated value(S)
        this.xics = deduplicatedXICs;
    }


    // TODO: Add function to serialize & unserialize here

    /**
     * Given a integer representing the window to be searched, deserialize the IsolationWindow.
     * @param serializedDataDir String representing the location of the serializedData
     * @param windowIndex int representing the window to be searched.
     * @return an unserialized IsolationWindow object.
     */
    public static IsolationWindow parseSerializedIsolationWindowFile(String serializedDataDir, int windowIndex){
        IsolationWindow deserializedIsolationWindowTrails = null;

        try {
            FileInputStream fileIn = new FileInputStream(serializedDataDir + windowIndex + "_window.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            deserializedIsolationWindowTrails = (IsolationWindow) in.readObject();
            in.close();
            fileIn.close();

        } catch (IOException | ClassNotFoundException $e){
            return null;
        }
        return deserializedIsolationWindowTrails;
    }


    public static void writeIsolationWindowPair(ArrayList<Pair<Double, Double>> isolationWindows, String pairOutfile) throws IOException {
        FileWriter fstream = new FileWriter(pairOutfile);
        BufferedWriter info = new BufferedWriter(fstream);
        int num_windows = isolationWindows.size();

        System.out.println("Writing IsolationWindow Ranges ...");

        for (int cur_iso = 0; cur_iso < num_windows; ++cur_iso) {
            info.write("START\n");
            info.write(
                    Double.toString(isolationWindows.get(cur_iso).getKey())
                            + " "
                            + Double.toString(isolationWindows.get(cur_iso).getValue())
                            + "\n");
//            info.write("END\n"); // TODO: refactor this if necessary
        }
        info.close();
        fstream.close();
    }

    public void writeIsolationWindowXICsToTSV(
            IsolationWindow isolationWindow, String workingDirectory, int windowIndex) throws IOException{
        FileWriter fileWriter = new FileWriter(workingDirectory + windowIndex + "_window.tsv");
        BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);

        // Write the header
        // maxRT, maxMZ, maxINT, rts, mzs, ints
        bufferedWriter.write("maxRT" + "\t" + "maxMZ" + "\t" + "maxINT" + "\t" + "rts" + "\t" + "mzs" + "\t" + "ints" + "\t" + "peakSum" + "\t" + "peakArea" + "\n");

        // For each XIC in XICs
        for (XIC xic : isolationWindow.xics) {
            // Write maxRT, maxMZ, maxINT
            bufferedWriter.write(Double.toString(xic.getRtAtMaxIntensity()) + "\t");
            bufferedWriter.write(Double.toString(xic.getMZAtMaxIntensity()) + "\t");
            bufferedWriter.write(Double.toString(xic.getMaxIntensity()) + "\t");

            // Write Retention Times
            StringBuilder stringBuilder = new StringBuilder();
            for (Double rts: xic.getRetentionTimes()){
                stringBuilder.append(Double.toString(rts)).append(",");
            }
            bufferedWriter.write(stringBuilder.toString() + "\t");

            // Write MZs
            stringBuilder = new StringBuilder();
            for (Double mzs: xic.getMassChargeRatios()){
                stringBuilder.append(Double.toString(mzs)).append(",");
            }
            bufferedWriter.write(stringBuilder.toString() + "\t");

            // Write INTs
            stringBuilder = new StringBuilder();
            for (Double ints: xic.getIntensities()){
                stringBuilder.append(Double.toString(ints)).append(",");
            }
            bufferedWriter.write(stringBuilder.toString() + "\t");

            // Write PEAKSUM
            bufferedWriter.write(Double.toString(xic.getPeakSum()) + "\t");

            // Write PEAKAREA
            bufferedWriter.write(Double.toString(xic.getPeakArea()) + "\n");
        }

        bufferedWriter.close();
        fileWriter.close();
    }


    /**
     * Given an ArrayList of XIC, subset the ArrayList by RT and some tolerance in minutes.
     * @return List of XIC objects that are extracted from the superset trails
     */
    public ArrayList subsetXICBasedOnRT(XIC baseXIC, double rtToleranceInMinutes){
        XIC upperRTXIC = new XIC(0d, baseXIC.getRtAtMaxIntensity() + rtToleranceInMinutes);
        XIC lowerRTXIC = new XIC(0d, baseXIC.getRtAtMaxIntensity() - rtToleranceInMinutes);

        /*
        int correctUpper = Collections.binarySearch(this.xics, upperRTXIC, Comparator.comparingDouble(XIC::getRtAtMaxIntensity));
        int correctLower = Collections.binarySearch(this.xics, lowerRTXIC, Comparator.comparingDouble(XIC::getRtAtMaxIntensity));
        correctUpper = parseInsertionPoint(correctUpper);
        correctLower = parseInsertionPoint(correctLower);
        */

        int upperIndex = Collections.binarySearch(this.xics, upperRTXIC, XIC::compareRtXIC);
        int lowerIndex = Collections.binarySearch(this.xics, lowerRTXIC, XIC::compareRtXIC);
        upperIndex = parseInsertionPoint(upperIndex);
        lowerIndex = parseInsertionPoint(lowerIndex);

        return new ArrayList<> (this.xics.subList(lowerIndex, upperIndex));
    }
    /**
     * Process the (possibly negative) index returned by Collections.binarySearch.
     * @param index int index of the binary searched item.
     * @return the index of the searched item OR the index to insert the new item.
     */
    private static int parseInsertionPoint(int index){
        if (index < 0){
            return (-1 * (index) - 1);
        } else
        {
            return index;
        }
    }

    public static ArrayList<edu.uw.waterlooms.entity.Pair> parseWindowRanges(String file) throws FileNotFoundException {
        try {
            ArrayList<edu.uw.waterlooms.entity.Pair> isolationWindows = new ArrayList<>();
            FileInputStream inputStream = new FileInputStream(file);
            Scanner sc = new Scanner(inputStream, "UTF-8");
            while (sc.hasNextLine()) {
                String line = sc.nextLine();
                if (line.contains("START")) {
                    line = sc.nextLine();
                    String[] range = line.split("\\s+");
                    double mzLow = Double.parseDouble(range[0]);
                    double mzHigh = Double.parseDouble(range[1]);
                    edu.uw.waterlooms.entity.Pair newWindow = new edu.uw.waterlooms.entity.Pair(mzLow, mzHigh);
                    isolationWindows.add(newWindow);
                }
            }
            return isolationWindows;
        } catch (FileNotFoundException $e) {

        }
        return new ArrayList<>();
    }
}


