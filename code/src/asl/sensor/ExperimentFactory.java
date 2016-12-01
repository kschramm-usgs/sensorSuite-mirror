package asl.sensor;

public class ExperimentFactory {
  
  // simple factory method for creating test data generators

  public static Experiment createExperiment(
      ExperimentEnum exp, 
      double[][] data){
    switch(exp){
    case ORTHO:
      return new OrthogonalExperiment(data);
    case NOISE:
      return new NoiseExperiment(data);
    case RGAIN:
      return new GainExperiment(data);
    default:
      throw new IllegalArgumentException("Invalid enum type specified");
    }
  }
  
}
