package asl.sensor;

/**
 * Enumerated type defining each kind of test, for GUI and test factory
 * @author akearns
 *
 */
public enum ExperimentEnum {

  // if adding a new test to this, make sure to also create a new extender for
  // experiment and a corresponding constructor call in experimentfactory
  
  NOISE ("Self-noise", false),
  RGAIN ("Relative gain", true),
  STCAL ("Step calibration", true),
  ORTHO ("Orthogonality", false);

  private String name;
  public boolean twoInput; 
    // not my favorite way of doing this but it's an easy modification to make
  
  ExperimentEnum(String name, boolean twoInput) {
    this.name = name;
    this.twoInput = twoInput;
  }
  
  public String getName() {
    return name;
  }
  
  public boolean isTwoInput() {
    return twoInput;
  }
  
}
