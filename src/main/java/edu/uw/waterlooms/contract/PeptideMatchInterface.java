package edu.uw.waterlooms.contract;

import edu.uw.waterlooms.entity.IsotopeFeature;
import edu.uw.waterlooms.match.Peptide;

import java.util.List;

public interface PeptideMatchInterface {

  // TODO: Add a tuple class here

  /**
   * Contract for outputting a set of Peptide-IsotopeFeature matches.
   *
   * @param peptides a list of peptides to be matched.
   * @param isotopeFeaturese a list of Isotope Features to be matched against.
   * @return a list of tuples containing Peptide-IsotopeFeature matches.
   */
  public void matchPeptidesToIsotopeFeatures(List<Peptide> peptides, List<IsotopeFeature> isotopeFeaturese);
}
