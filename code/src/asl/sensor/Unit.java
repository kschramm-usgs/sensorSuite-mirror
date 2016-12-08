package asl.sensor;

/**
 * Enum for specifiying unit types, as in values parsed from a response file
 * @author akearns
 *
 */
public enum Unit {

  /**
   * Specifies that a given piece of data has units of m/s
   */
  VELOCITY,
  /**
   * Specifies that a given piece of data has units of m/(s^2)
   */
  ACCELERATION;
  
}
