package edu.uw.waterlooms.entity;

import java.util.List;

public class IsotopeFeature {
  // TODO: Can a IsotopeFeature have >1 precursor?
  private XIC precursor;

  private List<XIC> trails;

  /**
   * Constructor for an IsotopeFeature.
   *
   * @param trails list of XIC/Trail(s) to add to this IsotopeFeature.
   */
  public IsotopeFeature(List<XIC> trails) {
    this.trails = trails;
  }

  /**
   * Getter for all the XIC/Trail(s) assigned to this IsotopeFeature.
   *
   * @return List<XIC>
   */
  public List<XIC> getTrails() {
    return this.trails;
  }

  /**
   * Add an XIC/Trail to the IsotopeFeature.
   *
   * @param trail XIC trail to be clustered to the IsotopeFeature.
   */
  public void addTrail(XIC trail) {
    this.trails.add(trail);
  }

  /**
   * Setter for IsotopeFeature.
   *
   * @param trails list of XIC/Trails to represent the IsotopeFeature.
   */
  public void setTrails(List<XIC> trails) {
    this.trails = trails;
  }
}
