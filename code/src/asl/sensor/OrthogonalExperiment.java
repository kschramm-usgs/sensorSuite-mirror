package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class OrthogonalExperiment extends Experiment {

  public OrthogonalExperiment(XYSeriesCollection dataIn) {
    super(dataIn);
    xAxisTitle = "Ortho X Axis (units)";
    yAxisTitle = "Ortho Y Axis (units)";
  }

  @Override
  XYSeriesCollection backend(XYSeriesCollection dataIn) {
    // TODO Auto-generated method stub
    return dataIn;
  }

}
