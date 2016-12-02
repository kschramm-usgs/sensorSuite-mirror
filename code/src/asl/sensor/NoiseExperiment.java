package asl.sensor;

import org.jfree.data.time.TimeSeriesCollection;

public class NoiseExperiment extends Experiment {

  public NoiseExperiment(TimeSeriesCollection dataIn) {
    super(dataIn);
    xAxisTitle = "Period (s)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
  }

  @Override
  TimeSeriesCollection backend(TimeSeriesCollection dataIn) {
    // TODO Auto-generated method stub
    return dataIn;
  }

}
