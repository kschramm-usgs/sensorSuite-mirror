package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class NoiseExperiment extends Experiment {

  public NoiseExperiment(XYSeriesCollection dataIn) {
    super(dataIn);
    xAxisTitle = "Period (s)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
  }

  @Override
  XYSeriesCollection backend(XYSeriesCollection dataIn) {
    // TODO Auto-generated method stub
    return dataIn;
  }

}
