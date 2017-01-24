package asl.sensor;

import java.awt.Font;
import java.text.SimpleDateFormat;
import java.util.List;
import java.util.TimeZone;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class StepExperiment extends Experiment {

  double f, h; //corner and damping of output (uncorrected)
  double fCorr, hCorr; // fit parameters to turn output into cal input
  
  public StepExperiment() {
    super();
    xAxisTitle = "Time (s)";
    yAxisTitle = "Counts (normalized)";
    xAxis = new DateAxis(xAxisTitle);
    SimpleDateFormat sdf = new SimpleDateFormat("HH:mm");
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    xAxis.setLabel("UTC Time");
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
  }
  
  public double[] getCornerAndDamping() {
    return new double[]{f, h};
  }
  
  @Override
  void backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    // assume that the first block is the raw step calibration
    // the raw calibration is defined as not having an associated response
    DataBlock stepCalRaw = ds.getXthLoadedBlock(1);
    // but we want the response and the data of the cal result
    int outIdx = ds.getXthFullyLoadedIndex(1);
    
    // if someone did load a raw cal with the response, then we wouldn't
    // get a different block with the second call above, so we get the 
    // next loaded block/response pair
    if ( ds.getBlock(outIdx).getName().equals( stepCalRaw.getName() ) ) {
      outIdx = ds.getXthFullyLoadedIndex(2);
    }
    

    DataBlock db = ds.getBlock(outIdx);
    long interval = db.getInterval();
    InstrumentResponse ir = ds.getResponse(outIdx);
    Complex pole = ir.getPoles().get(0);
    
    f = 1. / (2 * Math.PI / pole.abs() ); // corner frequency
    h = Math.abs( pole.getReal() / pole.abs() ); // damping
    Unit respUnit = ir.getUnits();
    List<Double> gain = ir.getGain();
    InstrumentResponse appxResp = new InstrumentResponse(f, h, respUnit, gain);
    
    // get FFT of datablock timeseries, apply response to input
    FFTResult fft = FFTResult.simpleFFT(db);
    
    Complex[] respValues = appxResp.applyResponseToInput( fft.getFreqs() );
    
    Complex maxVal = Complex.ZERO; // negative infinity
    for (Complex respVal : respValues) {
      if ( respVal.abs()  > maxVal.abs()  && respVal.abs() != Double.NaN) {
        maxVal = respVal;
      }
    }
    
    Complex[] fftValues = fft.getFFT();
    
    Complex[] correctedValues = new Complex[fftValues.length];
    // don't let denominator be zero
    System.out.println(fftValues[0]);
    System.out.println(respValues[0]);
    System.out.println(maxVal.multiply(0.008));
    correctedValues[0] = Complex.ONE;
    for (int i = 1; i < correctedValues.length; ++i) {
      Complex numer = fftValues[i].multiply( respValues[i].conjugate() );
      Complex denom = respValues[i].multiply( respValues[i].conjugate() );
      if ( denom.abs() == 0. ) {
        denom = denom.add( maxVal.abs() * 0.008 );
      }
      correctedValues[i] = numer.divide(denom);
    }
    
    double[] toPlot = FFTResult.inverseFFT(correctedValues);
    long now = db.getStartTime();
    
    XYSeries xys = new XYSeries( db.getName() );
    for (double point : toPlot) {
      xys.add(now/1000, point);
      now += interval;
    }
    
    // next we'll want to find the parameters to fit the plots
    // to the inputted data
    
    xySeriesData.addSeries(xys);
    
  }

}
