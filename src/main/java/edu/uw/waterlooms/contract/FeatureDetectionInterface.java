package edu.uw.waterlooms.contract;

import edu.uw.waterlooms.entity.IsotopeFeature;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.msutil.OpenMzxml;
import java.util.List;

public interface FeatureDetectionInterface {

  /**
   * Contract for outputting a list of XIC/Trails given a mzXML input.
   *
   * @param inFile the mzXML file to be parsed.
   * @return a list XIC/Trails.
   */
  public List<XIC> detectFeatures(OpenMzxml inFile);
}
