package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class GainExperiment extends Experiment {

  public GainExperiment(XYSeriesCollection dataIn) {
    super(dataIn);
    xAxisTitle = "Gain X Axis (units)";
    yAxisTitle = "Gain Y Axis (units)";
  }

  @Override
  XYSeriesCollection backend(XYSeriesCollection dataIn) {
    // TODO Auto-generated method stub
    return dataIn;
  }

}
