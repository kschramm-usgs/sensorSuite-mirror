package asl.sensor;

import java.awt.Font;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

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
    return xAxis;
  }
  
  public String getFreqTitle() {
    return xAxisTitle;
  }
  
  @Override
  void backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
 // TODO move to stepexperiment
    DataBlock db = ds.getBlock(0);
    long interval = db.getInterval();
    InstrumentResponse ir = ds.getResponse(0);
    Complex pole = ir.getPoles().get(0);
    // WE WILL WANT TO KEEP TRACK OF THESE VALUES
    double f = 1. / (2 * Math.PI / pole.abs() ); // corner frequency
    double h = Math.abs( pole.getReal() / pole.abs() ); // damping
    
    InstrumentResponse appxResp = new InstrumentResponse(f, h);
    
    // get FFT of datablock timeseries, apply response to input
    FFTResult fft = FFTResult.simpleFFT(db);
    
    Complex[] respValues = appxResp.applyResponseToInput( fft.getFreqs() );
    
    double maxVal = respValues[0].abs();
    for (Complex respVal : respValues) {
      if ( respVal.abs()  > maxVal ) {
        maxVal = respVal.abs();
      }
    }
    
    Complex[] fftValues = fft.getFFT();
    
    Complex[] correctedValues = new Complex[fftValues.length];
    // don't let denominator be zero
    for (int i = 0; i < correctedValues.length; ++i) {
      Complex numer = fftValues[i].multiply( respValues[i].conjugate() );
      Complex denom = respValues[i].multiply( respValues[i].conjugate() );
      denom = denom.add(0.008 * maxVal);
      correctedValues[i] = numer.divide(denom);
    }
    
    double[] toPlot = FFTResult.inverseFFT(correctedValues);
    long now = 0;
    
    XYSeries xys = new XYSeries( db.getName() );
    for (double point : toPlot) {
      xys.add(point, now);
      now += interval;
    }
    
    // next we'll want to find the parameters to fit the plots
    // to the inputted data
    
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(xys);
    
  }

}
