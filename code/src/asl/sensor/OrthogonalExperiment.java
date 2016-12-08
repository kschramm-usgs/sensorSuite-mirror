package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class OrthogonalExperiment extends Experiment {

  public OrthogonalExperiment() {
    super();
    xAxisTitle = "Ortho X Axis (units)";
    yAxisTitle = "Ortho Y Axis (units)";
  }

  @Override
  XYSeriesCollection backend(DataStore ds) {
    // TODO Auto-generated method stub
    return null;
  }

}
