package asl.sensor;

import java.awt.Font;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;

public class GainExperiment extends Experiment {

  private static int getPeakIndex(FFTResult fft) {
    
    Complex[] timeSeries = fft.getFFT();
    double[] freqs = fft.getFreqs();
    
    double max = Double.NEGATIVE_INFINITY;
    int index = 0;
    for (int i = 0; i < timeSeries.length; ++i) {
      Complex temp = timeSeries[i].multiply(Math.pow(2*Math.PI*freqs[i],4));
      double result = 10*Math.log10( temp.abs() );
      if ( result < Double.POSITIVE_INFINITY && result > max ) {
        max = result;
        index = i;
      }
    }
    return index;
  }
  
  private static double mean(FFTResult psd, int lower, int higher) {
    double result = 0;
    int range = higher-lower;
    
    Complex[] data = psd.getFFT();
    
    for (int i = lower; i <= higher; ++i) {
      if ( data[i].abs() >= Double.POSITIVE_INFINITY ) {
        continue;
      }
      result += data[i].abs();
    }
    
    return result/range; // since result is a double, no cast needed?
    
  }
  
  private static double sdev(
      FFTResult fft1, FFTResult fft2, 
      double meanRatio, 
      int lower, int higher) {
    Complex[] density1 = fft1.getFFT();
    Complex[] density2 = fft2.getFFT();
    double sigma = 0.;
    for (int i = lower; i <= higher; ++i) {
      double value1 = density1[i].abs();
      double value2 = density2[i].abs();
      
      if (value1 >= Double.POSITIVE_INFINITY || 
          value2 >= Double.POSITIVE_INFINITY) {
        continue;
      }
      
      sigma += Math.pow( (value1 / value2) - meanRatio, 2 );
    }
    
    return Math.sqrt(sigma);
  }
  
  private List<FFTResult> fftResults;
  
  private double ratio;

  private double sigma;
  
  public GainExperiment() {
    super();
    freqAxisTitle = "Period (s)";
    xAxisTitle = freqAxisTitle;
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
  void backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    fftResults = new ArrayList<FFTResult>( Arrays.asList(psd) );
    
    freqSpace = false;
    
    DataBlock[] dataIn = ds.getData();
    // used to get values related to range, etc.
    
    xySeriesData.setAutoWidth(true);
    
    addToPlot(dataIn, psd, freqSpace, xySeriesData, this);
    
    xySeriesData.addSeries( FFTResult.getLowNoiseModel(freqSpace, this) );
    
  }
  
  @Override
  public String[] getBoldSeriesNames() {
    return new String[]{"NLNM"};
  }

  private static int[] getRange(double[] freqs, double[] freqBoundaries) {
    
    double lowFrq = freqBoundaries[0];
    double highFrq = freqBoundaries[1];
    
    int[] indices = new int[2];
    
    for (int i = 1; i < freqs.length; ++i) {
      if (freqs[i] == lowFrq || 
          (freqs[i-1] < lowFrq && freqs[i] > lowFrq) ) {
        indices[0] = i;
      }
    }
    
    indices[1] = freqs.length-1;
    
    for (int i = indices[0]; i < freqs.length-1; ++i) {
      if (freqs[i] == highFrq || 
          (freqs[i] < highFrq && freqs[i+1] > highFrq) ) {
        indices[1] = i+1;
      }
    }
    
    return indices;
  }
  
  public double[] getStatsFromFreqs(int idx0, int idx1, 
      double lowFreq, double highFreq) {
    
    FFTResult plot0 = fftResults.get(idx0);
    
    double[] freqBoundaries = new double[2];
    freqBoundaries[0] = Math.min(lowFreq, highFreq);
    freqBoundaries[1] = Math.max(lowFreq, highFreq);
    
    int[] indices = getRange( plot0.getFreqs(), freqBoundaries );
    
    return getStatsFromIndices(idx0, idx1, indices[0], indices[1]);
  }
  
  private double[] getStatsFromIndices(int idx0, int idx1,
      int lowInd, int highInd) {
    
    FFTResult plot0 = fftResults.get(idx0);
    FFTResult plot1 = fftResults.get(idx1);
    
    double mean0 = mean(plot0, lowInd, highInd);
    // since both datasets must have matching interval, PSDs have same freqs
    double mean1 = mean(plot1, lowInd, highInd);
    
    // double MIN_VALUE field is effectively java's machine epsilon
    
    ratio = (mean0+Double.MIN_VALUE) / (mean1+Double.MIN_VALUE); 
      // prevent division by 0
    
    sigma = sdev(plot0, plot1, ratio, lowInd, highInd);
    
    // calculate ratio and sigma over the range
    
    return new double[]{ratio, sigma};
  }

}
