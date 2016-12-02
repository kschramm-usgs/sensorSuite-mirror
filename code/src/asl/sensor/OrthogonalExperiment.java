package asl.sensor;

import org.jfree.data.time.TimeSeriesCollection;

public class OrthogonalExperiment extends Experiment {

  public OrthogonalExperiment(TimeSeriesCollection dataIn) {
    super(dataIn);
    xAxisTitle = "Ortho X Axis (units)";
    yAxisTitle = "Ortho Y Axis (units)";
  }

  @Override
  TimeSeriesCollection backend(TimeSeriesCollection dataIn) {
    // TODO Auto-generated method stub
    return dataIn;
  }

}
