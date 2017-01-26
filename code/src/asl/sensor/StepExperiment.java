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
    xAxis = new NumberAxis(xAxisTitle);
    //SimpleDateFormat sdf = new SimpleDateFormat("HH:mm");
    //sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    //xAxis.setLabel("UTC Time");
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
    
    // resulting inverse fft should be this length
    int trimLength = stepCalRaw.getData().size();

    // get data of the result of the step calibration
    DataBlock sensorOutput = ds.getBlock(outIdx);
    long interval = sensorOutput.getInterval();
    InstrumentResponse ir = ds.getResponse(outIdx);
    Complex pole = ir.getPoles().get(0);
    
    f = 1. / (2 * Math.PI / pole.abs() ); // corner frequency
    h = Math.abs( pole.getReal() / pole.abs() ); // damping
    //Unit respUnit = ir.getUnits();
    //List<Double> gain = ir.getGain();
    //InstrumentResponse apxResp = new InstrumentResponse(f, h, respUnit, gain);
    
    // get FFT of datablock timeseries, apply response to input
    FFTResult fft = FFTResult.simpleFFT(stepCalRaw);
    
    // term inside the square root in the calculations of p1, p2
    // (h^2-1)
    Complex tempResult = new Complex( Math.pow(h, 2) - 1 );
    
    double omega = 2 * Math.PI * f; // omega_0
    

    // casting h as a complex number to enable these calculations
    Complex hCast = new Complex(h);
    
    // - (h + sqrt(h^2-1))
    Complex pole1 = hCast.add( tempResult.sqrt() ).multiply(-1);
    pole1.multiply(omega);
    
    // - (h - sqrt(h^2-1))
    Complex pole2 = hCast.subtract( tempResult.sqrt() ).multiply(-1);
    pole2.multiply(omega);
    
    Complex[] fftValues = fft.getFFT();
    double[] freqs = fft.getFreqs();
    
    Complex[] respPerFreq = new Complex[fftValues.length]; // array of resps
    double max = 0.0;
    // don't let denominator be zero
    respPerFreq[0] = Complex.ONE;
    for (int i = 1; i < respPerFreq.length; ++i) {
      
      // 2*pi*i*f
      Complex factor = new Complex(0, 2*Math.PI*freqs[i]); 
      
      // (2*pi*i*f - p1) * (2*pi*f*i - p2)
      Complex denom = factor.subtract(pole1).multiply( factor.subtract(pole2) );
      
      Complex resp = factor.divide(denom);
      
      // RESP = fft(x)*k / ((k-p1)(k-p2))
      // where k = 2*pi*i*f, x is the input signal, f is FFT frequencies
      respPerFreq[i] = 
          fftValues[i].multiply(resp);
      
      if (respPerFreq[i].abs() > max) {
        max = respPerFreq[i].abs();
      }
      
    }
    
    Complex[] stepFFT = FFTResult.simpleFFT(sensorOutput).getFFT();
    Complex[] correctedValues = new Complex[stepFFT.length];
    
    for (int i = 0; i < stepFFT.length; ++i) {
      // the conjugate of the response, used twice in deconvolve function
      Complex conjResp = respPerFreq[i].conjugate();
      double aboveZero = 0.008*max;  // term to keep the denominator above 0
      
      // resp * conj(resp) + 0.008(max(|resp|))
      Complex denom = respPerFreq[i].multiply(conjResp).add(aboveZero);
      
      // deconvolving the response
      // fft * conj(resp) / (resp * conjResp)+0.008(max(|resp|))
      correctedValues[i] = stepFFT[i].multiply(conjResp).divide(denom);
    }
    
    double[] toPlot = FFTResult.inverseFFT(correctedValues, trimLength);
    
    /*
    toPlot = new double[stepFFT.length];
    for (int i = 0; i < correctedValues.length; ++i) {
      toPlot[i] = correctedValues[i].getReal();
    }
    */
    
    long start = 0L; // was db.startTime();
    long now = start;
    XYSeries xys = new XYSeries("FFT(STEP)/RESP");
    XYSeries scs = new XYSeries( stepCalRaw.getName() ); 
    for (double point : toPlot) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      xys.add(seconds, point);
      now += interval;
    }
    now = start;
    for (Number point : stepCalRaw.getData() ) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      scs.add(seconds, point);
      now += interval;
    }
    
    // next we'll want to find the parameters to fit the plots
    // to the inputted data
    // xySeriesData.addSeries(scs);
    xySeriesData.addSeries(xys);
    
  }

}
