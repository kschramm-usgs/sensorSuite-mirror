package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class GainExperiment extends Experiment {

  public GainExperiment() {
    super();
    xAxisTitle = "Gain X Axis (units)";
    yAxisTitle = "Gain Y Axis (units)";
  }

  @Override
  XYSeriesCollection backend(DataStore ds) {
    // TODO Auto-generated method stub
    return null;
  }

}
