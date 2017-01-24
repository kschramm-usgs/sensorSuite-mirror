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
    //Unit respUnit = ir.getUnits();
    //List<Double> gain = ir.getGain();
    //InstrumentResponse apxResp = new InstrumentResponse(f, h, respUnit, gain);
    
    // get FFT of datablock timeseries, apply response to input
    FFTResult fft = FFTResult.simpleFFT(db);
    
    // term inside the square root in the calculations of p1, p2
    Complex tempResult = new Complex( Math.pow(h, 2) - 1 );
    
    double omega = 2 * Math.PI * f; // omega_0
    
    Complex pole1 = tempResult.sqrt().add(h).multiply(-1);
    pole1.multiply(omega);
    
    Complex pole2 = tempResult.sqrt().subtract(h).multiply(-1);
    pole2.multiply(omega);
    
    /*
    Complex[] respValues = apxResp.applyResponseToInput( fft.getFreqs() );
    
    Complex maxVal = Complex.ZERO; // negative infinity
    for (Complex respVal : respValues) {
      if ( respVal.abs()  > maxVal.abs()  && respVal.abs() != Double.NaN) {
        maxVal = respVal;
      }
    }
    */
    
    Complex[] fftValues = fft.getFFT();
    double[] freqs = fft.getFreqs();
    
    Complex[] correctedValues = new Complex[fftValues.length];
    // don't let denominator be zero
    correctedValues[0] = Complex.ONE;
    for (int i = 1; i < correctedValues.length; ++i) {
      Complex factor = new Complex(0, 2*Math.PI*freqs[i]); // 2*pi*i*f
      
      // fft*(2*pi*i*f-0)
      Complex numer = fftValues[i].multiply(factor);

      // (2*pi*i*f-p1)*(2*pi*i*f-p2)
      Complex denom = factor.subtract(pole1);
      denom = denom.multiply( factor.subtract(pole2) );
      
      correctedValues[i] = numer.divide(denom);
    }
    
    double[] toPlot = FFTResult.inverseFFT(correctedValues);
    long start = db.getStartTime();
    System.out.println("start time: " + start);
    long now = start;
    XYSeries xys = new XYSeries( db.getName() );
    XYSeries scs = new XYSeries( stepCalRaw.getName() ); 
    for (double point : toPlot) {
      xys.add(now/1000, point);
      now += interval;
    }
    now = start;
    for (Number point : stepCalRaw.getData() ) {
      scs.add(now/1000, point);
      now += interval;
    }
    
    // next we'll want to find the parameters to fit the plots
    // to the inputted data
    // xySeriesData.addSeries(scs);
    xySeriesData.addSeries(xys);
    
  }

}
