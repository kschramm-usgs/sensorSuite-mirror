package asl.sensor;

public class ExperimentPanelFactory {

  public static ExperimentPanel createPanel(
      ExperimentEnum exp){
    switch(exp){
    case ORTHO:
      return new OrthogonalPanel(exp);
    case NOISE:
      return new NoisePanel(exp);
    case RGAIN:
      return new GainPanel(exp);
    default:
      // this shouldn't happen unless someone added to the enum
      // and forgot to follow-through on implementation
      throw new IllegalArgumentException("Invalid enum type specified");
    }
  }
  
}
