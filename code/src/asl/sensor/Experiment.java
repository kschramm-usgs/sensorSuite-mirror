package asl.sensor;

/**
 * Enumerated type defining each kind of test, for GUI and test factory
 * @author akearns
 *
 */
public enum Experiment {

  ORTHO ("Orthogonality"),
  RGAIN ("Relative gain"),
  NOISE ("Self-noise");
  
  private String name;
  
  Experiment(String name) {
    this.name = name;
  }
  
  public String getName() {
    return name;
  }
  
}
