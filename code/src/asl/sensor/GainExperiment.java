package asl.sensor;

import org.jfree.data.time.TimeSeriesCollection;

public class GainExperiment extends Experiment {

  public GainExperiment(TimeSeriesCollection dataIn) {
    super(dataIn);
  }

  @Override
  TimeSeriesCollection backend(TimeSeriesCollection dataIn) {
    // TODO Auto-generated method stub
    return dataIn;
  }

}
