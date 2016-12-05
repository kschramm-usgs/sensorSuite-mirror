package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class ExperimentFactory {
  
  // simple factory method for creating test data generators

  public static Experiment createExperiment(
      ExperimentEnum exp, 
      XYSeriesCollection data){
    switch(exp){
    case ORTHO:
      return new OrthogonalExperiment(data);
    case NOISE:
      return new NoiseExperiment(data);
    case RGAIN:
      return new GainExperiment(data);
    default:
      // this shouldn't happen unless someone added to the enum
      // and forgot to follow-through on implementation
      throw new IllegalArgumentException("Invalid enum type specified");
    }
  }
  
}
