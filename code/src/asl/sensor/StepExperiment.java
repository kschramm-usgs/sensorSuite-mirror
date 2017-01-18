package asl.sensor;

import java.awt.Font;

import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;

public class StepExperiment extends Experiment {

  public StepExperiment() {
    super();
    xAxisTitle = "Time (s)";
    yAxisTitle = "Counts (normalized)";
    xAxis = new NumberAxis(xAxisTitle);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRangeIncludesZero(false);
    yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
  }
  
  @Override
  public ValueAxis getFreqAxis() {
    throw new RuntimeException("X-axis is a timescale, not freq/period");
  }
  
  public String getFreqTitle() {
    throw new RuntimeException("X-axis is a timescale, not freq/period");
  }
  
  @Override
  void backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    // we'll need to take FFTs, but not the PSD
    // TODO Auto-generated method stub

    // what do we do here? get corner freq and damping from first pole
    
    // then take the FFT of the inputted data and multiply by response
    // (i.e., convolving time-series data with response function in time-dmn)
    // then de-convolve with the calculated response (from corner & damping)
    
    // after that do Levenberg-Marquardt to match to the input timeseries
    // giving us known parameters from response and calculated ones
    // calculated ones effectively matching the input
    
    
    
  }

}
