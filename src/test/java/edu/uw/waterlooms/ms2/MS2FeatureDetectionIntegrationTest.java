package edu.uw.waterlooms.ms2;

import edu.uw.waterlooms.entity.IsolationWindow;
import org.junit.Test;

import java.io.IOException;
import java.nio.file.FileSystems;

public class MS2FeatureDetectionIntegrationTest {

    // @Test
    public void ms2FeatureDetectionShouldReturnXICsForAllIsolationWindows(){
        // Working class to ensure no duplication(s)
        /* Given */
        String dataDirectory = FileSystems.getDefault().getPath("").toAbsolutePath().toString() + "/data/";
        String mzxmlFile = dataDirectory + "dia_data.mzXML";

        IsolationWindow window = null;

        MS2FeatureDetection ms2FeatureDetection = new MS2FeatureDetection(20, dataDirectory, mzxmlFile);
        try {
            window = ms2FeatureDetection.extractMS2XICForASingleIsolationWindow(0, true);
        } catch (IOException $exception){
            System.err.println($exception.getMessage());
        }

        assert window != null;
        window.removeDuplicateTrails();

        int removethis = 0;





        /* When */


        /* Then */

    }
}
