package asl.sensor;

import java.awt.Font;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;

public class GainExperiment extends Experiment {

  private double[] gainStage1;
  
  /**
   * Finds the maximum value of PSD plot curve
   * @param fft PSD calculation, including both FFT function and matching
   * frequencies
   * @return The index of the peak location
   */
  private int getPeakIndex(int idx) {
    
    FFTResult fft = fftResults.get(idx);
    
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
  
  public double[] getOctaveCenteredAtPeak(int idx) {
    
    int center = getPeakIndex(idx);
    double[] freqs = fftResults.get(idx).getFreqs();
    double peakFreq = freqs[center];
    
    double lowFreq = peakFreq / Math.sqrt(2);
    double highFreq = peakFreq * Math.sqrt(2);
    
    return new double[]{lowFreq, highFreq};
  }
  
  public double[] getStatsFromPeak(int idx0, int idx1) {
    double[] freqBounds = getOctaveCenteredAtPeak(idx0);
    double freq1 = freqBounds[0];
    double freq2 = freqBounds[1];
    return getStatsFromFreqs(idx0, idx1, freq1, freq2);
  }
  
  /**
   * Gets the indices denoting the inclusive range of a frequency bound on
   * a list of input frequencies
   * @param freqs Frequency series resulting from FFT calculation
   * (must be pre-sorted)
   * @param freqBoundaries Array of size 2 denoting frequency 
   * upper and lower bounds (in Hz), lower bound first
   * @return Array of size 2 denoting the lower and upper indices of the 
   * sub-range of frequencies including the given boundaries
   */
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
  
  /**
   * Get the mean of the PSD calculation within the specified range
   * @param psd PSD calculation (frequency space / FFT)
   * @param lower Starting index of window
   * @param higher Ending index of window
   * @return Mean of datapoints in the specific range
   */
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
  
  /**
   * PSD Standard deviation calculation given mean ratio and specified range
   * (Get the result for mean calculation using the mean function before
   * calling this calculation; does not call the function directly)
   * @param fft1 first PSD (numerator) of relative gain standard deviation
   * @param fft2 second PSD (denominator)
   * @param meanRatio mean(fft1)/mean(fft2)
   * @param lower Starting index of window (should be same used in mean calc)
   * @param higher Ending index of window (should be same used in mean calc)
   * @return double corresponding to standard deviation ratio over the window
   */
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
  
  /**
   * Constructor for the gain experiment (setting axes)
   * Currently does not include setting to set x-axis to display in freq. space
   */
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
  
  /**
   * Populate XYSeriesCollection with all input data (will be trimmed on plot)
   */
  @Override
  void backend(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    gainStage1 = new double[DataStore.FILE_COUNT];
    InstrumentResponse[] resps = ds.getResponses();
    for (int i = 0; i < gainStage1.length; ++i) {
      if (resps[i] == null) {
        continue;
      }
      double[] gains = resps[i].getGain();
      gainStage1[i] = gains[1];
    }
    
    fftResults = new ArrayList<FFTResult>( Arrays.asList(psd) );
    
    freqSpace = false;
    
    DataBlock[] dataIn = ds.getData();
    // used to get values related to range, etc.
    
    xySeriesData.setAutoWidth(true);
    
    addToPlot(dataIn, psd, freqSpace, xySeriesData, this);
    
    xySeriesData.addSeries( FFTResult.getLowNoiseModel(freqSpace, this) );
    
  }

  @Override
  /**
   * Specifies plotting the NLNM curve in bold
   */
  public String[] getBoldSeriesNames() {
    return new String[]{"NLNM"};
  }
  
  /**
   * Given indices to specific PSD data sets and frequency boundaries, gets
   * the mean and standard deviation ratios 
   * @param idx0 Index of numerator PSD
   * @param idx1 Index of denominator PSD
   * @param lowFreq Lower-bound of frequency window of PSD
   * @param highFreq Upper-bound of frequency window of PSD
   * @return 2-entry array of form {mean, standard deviation}
   */
  public double[] getStatsFromFreqs(int idx0, int idx1, 
      double lowFreq, double highFreq) {
    
    // TODO: check that indices are within a valid range
    
    FFTResult plot0 = fftResults.get(idx0);
    
    double[] freqBoundaries = new double[2];
    freqBoundaries[0] = Math.min(lowFreq, highFreq);
    freqBoundaries[1] = Math.max(lowFreq, highFreq);
    
    int[] indices = getRange( plot0.getFreqs(), freqBoundaries );
    
    return getStatsFromIndices(idx0, idx1, indices[0], indices[1]);
  }
  
  /**
   * Given indices to specific PSD data sets and indices to the corresponding
   * frequency boundaries, gets the mean and standard deviation ratios
   * @param idx0 Index of numerator PSD
   * @param idx1 Index of denominator PSD
   * @param lowInd Lower-bound index of PSDs' frequency array
   * @param highInd Upper-bound index of PSDs' frequency array
   * @return
   */
  private double[] getStatsFromIndices(int idx0, int idx1,
      int lowInd, int highInd) {
    
    // make sure lowInd really is the lower index
    if (lowInd > highInd) {
      int temp = highInd;
      highInd = lowInd;
      lowInd = temp;
    }
    
    FFTResult plot0 = fftResults.get(idx0);
    FFTResult plot1 = fftResults.get(idx1);
    
    double mean0 = mean(plot0, lowInd, highInd);
    // since both datasets must have matching interval, PSDs have same freqs
    double mean1 = mean(plot1, lowInd, highInd);
    
    // double MIN_VALUE field is effectively java's machine epsilon
    
    ratio = (mean0+Double.MIN_VALUE) / (mean1+Double.MIN_VALUE); 
      // prevent division by 0
    
    sigma = sdev(plot0, plot1, ratio, lowInd, highInd);
    
    double gain1 = gainStage1[idx0];
    double gain2 = gainStage1[idx1]/Math.sqrt(ratio);
    
    // calculate ratio and sigma over the range
    
    return new double[]{ratio, sigma, gain1, gain2};
  }

}
