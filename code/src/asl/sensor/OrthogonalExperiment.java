package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public class OrthogonalExperiment extends Experiment {

  public OrthogonalExperiment() {
    super();
    xAxisTitle = "Ortho X Axis (units)";
    yAxisTitle = "Ortho Y Axis (units)";
  }

  @Override
  XYSeriesCollection backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    // START OF COPY
    DataBlock[] dataIn = ds.getData();
    // InstrumentResponse[] responses = ds.getResponses();
    
    XYSeriesCollection plottable = new XYSeriesCollection();
    plottable.setAutoWidth(true);
    
    ExperimentPanel.addToPlot(dataIn, psd, plottable, freqSpace);
    
    return plottable;
    
  }

}
