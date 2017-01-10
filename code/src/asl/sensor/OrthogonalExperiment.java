package asl.sensor;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class OrthogonalExperiment extends Experiment {

  public OrthogonalExperiment() {
    super();
    xAxisTitle = "Ortho X Axis (units)";
    yAxisTitle = "Ortho Y Axis (units)";
  }

  @Override
  void backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    // TODO auto generated method stub
  }



}
