package asl.sensor;

import java.awt.Font;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;

public class GainExperiment extends Experiment {

  private static int getPeakIndex(Complex[] timeSeries, double[] freqs) {
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
      sigma += Math.pow( (value1 / value2) - meanRatio, 2 );
    }
    
    return Math.sqrt(sigma);
  }
  
  private FFTResult plot0;
  private FFTResult plot1;
  
  private double[] freqBoundaries = new double[2];
  
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
    
    plottingIndices = new int[]{0, 1}; // defaults
  }
  
  @Override
  void backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    freqSpace = false;
    
    DataBlock[] dataIn = ds.getData();
    // used to get values related to range, etc.
    plot0 = psd[ plottingIndices[0] ];
    plot1 = psd[ plottingIndices[1] ];
    
    getCenteredOctave();
    // sets freqBoundaries as a side effect
    
    xySeriesData.setAutoWidth(true);
    
    addToPlot(dataIn, psd, freqSpace, plottingIndices, xySeriesData, this);
    
    xySeriesData.addSeries( FFTResult.getLowNoiseModel(freqSpace, this) );
    
  }
  
  @Override
  public String[] getBoldSeriesNames() {
    return new String[]{"NLNM"};
  }
  
  private double[] getCenteredOctave() {
    
    // find the octave centered around frequency x, where peak of PSD is
    // that is, the limits are halfway between (x/2, x) and (x, 2x)
    // thus at the locations of 3x/4 and 3x/2
    
    int peak = getPeakIndex( plot0.getFFT(), plot0.getFreqs() );
    
    double[] freqs = plot0.getFreqs();
    
    double centerFrq = freqs[peak];
    
    freqBoundaries[0] = centerFrq * 3./4.; // 
    freqBoundaries[1] = centerFrq * 3./2.;
    
    int[] range = getRange(freqs, freqBoundaries);
    
    return getStatsFromIndices(range[0], range[1]);
    
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
  
  public double[] getStats() {
    return new double[]{ratio, sigma};
  }
  
  public double[] getFreqBoundaries() {
    return freqBoundaries;
  }
  
  public double[] getStatsFromFreqs(double lowFreq, double highFreq) {
    freqBoundaries[0] = Double.min(lowFreq, highFreq);
    freqBoundaries[1] = Double.max(lowFreq, highFreq);
    
    int[] indices = getRange( plot0.getFreqs(), freqBoundaries );
    
    return getStatsFromIndices(indices[0], indices[1]);
  }
  
  private double[] getStatsFromIndices(int lowInd, int highInd) {
    double mean0 = mean(plot0, lowInd, highInd);
    // since both datasets must have matching interval, PSDs have same freqs
    double mean1 = mean(plot1, lowInd, highInd);
    
    // double MIN_VALUE field is effectively java's machine epsilon
    
    ratio = (mean0+Double.MIN_VALUE) / (mean1+Double.MIN_VALUE); 
      // prevent division by 0
    
    sigma = sdev(plot0, plot1, ratio, lowInd, highInd);
    
    // calculate ratio and sigma over the range
    
    return getStats();
  }
  
  @Override
  public void setPlottingIndices(int[] plotBlocks) {
    if (plotBlocks.length != 2) {
      throw new RuntimeException("Relative gain requires exactly 2 plots");
    }
    plottingIndices = plotBlocks;
  }

}
