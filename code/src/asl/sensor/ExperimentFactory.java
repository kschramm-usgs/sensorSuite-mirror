package asl.sensor;

/**
 * Simple factory method for specifying a type of sensor test to calculate
 * @author akearns
 *
 */
public class ExperimentFactory {
  
  /**
   * Instantiates an Experiment from which the ExperimentPanel generates data
   * based on the type of Experiment given by the enum
   * @param exp Type of experiment (see the enum for details on kinds)
   * @return An implementation of the Experiment interface
   */
  public static Experiment createExperiment(
      ExperimentEnum exp){
    switch(exp){
    case ORTHO:
      return new OrthogonalExperiment();
    case NOISE:
      return new NoiseExperiment();
    case RGAIN:
      return new GainExperiment();
    case STCAL:
      return new StepExperiment();
    default:
      // this shouldn't happen unless someone added to the enum
      // and forgot to follow-through on implementation
      throw new IllegalArgumentException("Invalid enum type specified");
    }
  }
  
}
