package asl.sensor;

/**
 * Enumerated type defining each kind of test, for GUI and test factory
 * @author akearns
 *
 */
public enum ExperimentEnum {

  // if adding a new test to this, make sure to also create a new extender for
  // experiment and a corresponding constructor call in experimentfactory
  
  NOISE ("Self-noise", 3, 3),
  RGAIN ("Relative gain", 2, 2),
  STCAL ("Step calibration", 2, 1),
  ORTHO ("Orthogonality", 3, 3);

  private String name;
  private int blocksNeeded;
  private int fullDataNeeded;
    // not my favorite way of doing this but it's an easy modification to make
  
  ExperimentEnum(String name, int blocksNeeded, int fullDataNeeded) {
    this.name = name;
    this.blocksNeeded = blocksNeeded;
    this.fullDataNeeded = fullDataNeeded;
  }
  
  public int blocksNeeded() {
    return blocksNeeded;
  }
  
  public int fullDataNeeded() {
    return fullDataNeeded;
  }
  
  public String getName() {
    return name;
  }
  
}
