package asl.sensor;

import java.awt.Font;

import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeriesCollection;

public class GainExperiment extends Experiment {

  public GainExperiment() {
    super();
    xAxisTitle = "Period (s)";
    freqAxisTitle = "Frequency (Hz)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    xAxis = new LogarithmicAxis(xAxisTitle);
    freqAxis = new LogarithmicAxis(freqAxisTitle);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRangeIncludesZero(false);
    yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
  }
  
  @Override
  public String[] getBoldSeriesNames() {
    return new String[]{"NLNM"};
  }

  @Override
  XYSeriesCollection backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    // TODO Auto-generated method stub
    return null;
  }

}
