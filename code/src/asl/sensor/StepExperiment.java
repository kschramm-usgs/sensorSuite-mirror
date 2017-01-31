package asl.sensor;

import java.awt.Font;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;

public class StepExperiment extends Experiment {

  double f, h; //corner and damping of output (uncorrected)
  double fCorr, hCorr; // fit parameters to turn output into cal input
  
  public StepExperiment() {
    super();
    xAxisTitle = "Time (s)";
    yAxisTitle = "Counts (normalized)";
    xAxis = new NumberAxis(xAxisTitle);
    xAxis.setAutoRange(true);
    //SimpleDateFormat sdf = new SimpleDateFormat("HH:mm");
    //sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    //xAxis.setLabel("UTC Time");
    yAxis = new NumberAxis(yAxisTitle);
    // yAxis.setAutoRange(true);
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
    
    // use these terms to calculate the poles
    
    // term inside the square root in the calculations of p1, p2
    // (h^2-1)
    Complex tempResult = new Complex( Math.pow(h, 2) - 1 );
    
    double omega = 2 * Math.PI * f; // omega_0
    
    // casting h as a complex number to enable these calculations
    Complex hCast = new Complex(h);
    
    // - (h + sqrt(h^2-1))
    Complex pole1 = hCast.add( tempResult.sqrt() ).multiply(-1);
    pole1 = pole1.multiply(omega); 
      // cannot modify complex numbers "in-place"; they MUST be (re)assigned
    
    // - (h - sqrt(h^2-1))
    Complex pole2 = hCast.subtract( tempResult.sqrt() ).multiply(-1);
    pole2 = pole2.multiply(omega);
    
    // get FFT of datablock timeseries, deconvolve with response
    FFTResult sensorsFFT = FFTResult.simpleFFT(sensorOutput);
    Complex[] fftValues = sensorsFFT.getFFT();
    double[] freqs = sensorsFFT.getFreqs();
    
    // calculate the FFT of the response
    Complex[] respFFT = new Complex[fftValues.length]; // array of resps
    double max = 0.0;
    // don't let denominator be zero
    respFFT[0] = Complex.ONE;
    for (int i = 1; i < respFFT.length; ++i) {

      // replaced freqs[i] with 
      // 2*pi*i*f
      Complex factor = new Complex(0, 2*Math.PI*freqs[i]); 
      
      // (2*pi*i*f - p1) * (2*pi*f*i - p2)
      Complex denom = factor.subtract(pole1).multiply( factor.subtract(pole2) );
      
      respFFT[i] = factor.divide(denom);
      
      
      if (respFFT[i].abs() > max) {
        max = respFFT[i].abs();
      }
      
    }
    
    XYSeries freqTestPlot = new XYSeries("Freq plot");
    XYSeries respTestPlot = new XYSeries("Resp-applied plot");
    XYSeries fftResPlot = new XYSeries("Forward FFT result (absval)");
    for (int i = 2000; i <= freqs.length/2+10000; ++i) {
      freqTestPlot.add(i, freqs[i]);
      respTestPlot.add(i, respFFT[i].abs());
      fftResPlot.add(i, fftValues[i].abs());
    }
    
    Complex[] correctedValues = new Complex[fftValues.length];
    
    // deconvolving response is dividing fft(signal) by fft(response)
    
    for (int i = 0; i < fftValues.length; ++i) {
      // the conjugate of the response, used twice in deconvolution
      Complex conjResp = respFFT[i].conjugate();
      double aboveZero = 0.008*max;  // term to keep the denominator above 0
      
      // resp * conj(resp) + 0.008(max(|resp|))
      Complex denom = respFFT[i].multiply(conjResp).add(aboveZero);
      
      // deconvolving the response from output
      // fft * conj(resp) / (resp * conjResp)+0.008(max(|resp|))
      correctedValues[i] = fftValues[i].multiply(conjResp).divide(denom);
    }
    
    double[] toPlot = FFTResult.inverseFFT(correctedValues, trimLength);
    
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
    xySeriesData.addSeries(freqTestPlot);
    xySeriesData.addSeries(respTestPlot);
    xySeriesData.addSeries(fftResPlot);
    //xySeriesData.addSeries(scs);
    //xySeriesData.addSeries(xys);
    
  }

}
