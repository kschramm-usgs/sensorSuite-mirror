package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class NoiseExperiment extends Experiment {

  public NoiseExperiment() {
    super();
    xAxisTitle = "Period (s)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
  }

  @Override
  XYSeriesCollection backend(DataBlock[] dataIn) {
    // TODO Auto-generated method stub
    return null;
  }
  
  private void calcCrossPower(){
    
  }

}
