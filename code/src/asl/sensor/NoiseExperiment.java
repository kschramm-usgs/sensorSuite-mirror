package asl.sensor;

import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Produces the data for a self-noise test. Calculates PSD to get cross-power.
 * Based on code in the seedscan timeseries package, see
 * https://github.com/usgs/seedscan/tree/master/src/main/java/asl/timeseries
 * @author akearns, jholland 
 *
 */
public class NoiseExperiment extends Experiment {

  /**
   * Specifies the width of the cosine taper function used in windowing
   */
  private static final double TAPER_WIDTH = 0.10;
  
  /**
   * Instantiates a noise experiment -- axis titles and scales
   */
  public NoiseExperiment() {
    super();
    xAxisTitle = "Period (s)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    xAxis = new LogarithmicAxis(xAxisTitle); 
    yAxis = new NumberAxis(yAxisTitle);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
  }
  
  @Override
  public String[] getBoldSeriesNames() {
    return new String[]{"NLNM"};
  }
  /**
   * Generates power spectral density of each inputted file, and calculates
   * self-noise based on that result.
   * The overhead view is as follows: 
   * Take a window of size 1/4 incrementing through 1/16 of the data and
   * calculate the FFT of that region. Average these results together.
   * Apply the magnitude of the frequency response (relative to the FFT indices)
   * to that result and then take the complex conjugate. 
   * This produces the PSD plots.
   * Then, take the cross-powers of each of the terms (same calculation, but
   * multiply one result by the complex conjugate of the other), producing the
   * remaining terms for the formula for the self-noise results.
   */
  @Override
  XYSeriesCollection backend(DataStore ds) {
    
    DataBlock[] dataIn = ds.getData();
    InstrumentResponse[] responses = ds.getResponses();
    
    XYSeriesCollection plottable = new XYSeriesCollection();
    plottable.setAutoWidth(true);
    
    // TODO: make sure (i.e., when reading in data) that series lengths' match
    // rather than throwing the exceptions here
    Complex[][] spectra = new Complex[dataIn.length][];
    Complex[][] freqRespd = new Complex[dataIn.length][];
    double[] freqs = new double[1]; // initialize to prevent later errors
    
    long interval = dataIn[0].getInterval();
    int length = dataIn[0].getData().size();
    for (DataBlock data : dataIn) {
      if ( data.getInterval() != interval ) {
        throw new RuntimeException("Interval mismatch on datasets.");
      }
      if( data.getData().size() != length ) {
        // TODO: truncate data on read-in so this doesn't happen
        throw new RuntimeException("Length mismatch on datasets.");
      }
    }
    
    for (int i = 0; i < dataIn.length; ++i) {
      // don't calculate if all the data isn't in yet
      if( dataIn[i] == null ||  
          dataIn[i].getData().size() == 0 ||
          responses[i] == null ) {
        return new XYSeriesCollection();
        // we can't plot without all the data (certainly need responses loaded)
      }
    }
    
    
    // TODO: send to an external method so it's not duplicated?
    for (int i = 0; i < dataIn.length; ++i) {
      
      DataBlock data = dataIn[i];
      InstrumentResponse ir = responses[i];
      XYSeries powerSeries = new XYSeries("PSD "+data.getName() );
      
      FFTStruct fft = crossPower(data, data, ir, ir);
      
      Complex[] resultPSD = fft.getFFT();
      
      spectra[i] = resultPSD; // array of pii values
      freqs = fft.getFreqs();
      
      // TODO: find a way to refactor since this is technically redundant
      freqRespd[i] = ir.applyResponseToInput(freqs);
      
      for (int j = 0; j < freqs.length; ++j) {
        if (1/freqs[j] > 1.0E3) {
          continue;
        }
         
        powerSeries.add( 1/freqs[j], 10*Math.log10( resultPSD[j].abs() ) );
      }
     
      plottable.addSeries(powerSeries);
      
    }
    
    // spectra[i] is crosspower pii, now to get pij terms for i!=j
    FFTStruct fft = 
        crossPower(dataIn[0], dataIn[2], responses[0], responses[2]);
    Complex[] c13 = fft.getFFT();
    fft = crossPower(dataIn[1], dataIn[0], responses[1], responses[0]);
    Complex[] c21 = fft.getFFT();
    fft = crossPower(dataIn[1], dataIn[2], responses[1], responses[2]);
    Complex[] c23 = fft.getFFT();
    
    // WIP: use PSD results to get noise at each point see spectra
    XYSeries[] noiseSeriesArr = new XYSeries[dataIn.length];
    for(int j = 0; j < dataIn.length; ++j) {
      // initialize each xySeries with proper name for the data
      noiseSeriesArr[j] = new XYSeries( "Noise " + dataIn[j].getName() );
    }
    
    for (int i = 1; i < freqs.length; ++i) {
        if (1/freqs[i] > 1.0E3){
          continue;
        }
        Complex f1 = freqRespd[0][i];
        Complex f2 = freqRespd[1][i];
        Complex f3 = freqRespd[2][i];
        
        Complex p11 = spectra[0][i];
        Complex p22 = spectra[1][i];
        Complex p33 = spectra[2][i];

        Complex p13 = c13[i];
        Complex p21 = c21[i];
        Complex p23 = c23[i];
        
        // nii = pii - pij*hij
        Complex n11 = 
            p11.subtract(
                p21.multiply(p13).divide(p23) )
              .divide(f1);
        n11 = n11.multiply(Math.pow(2*Math.PI*freqs[i],2));
        
        Complex n22 = 
            p22.subtract(
                ( p23.conjugate() ).multiply(p21).divide( p13.conjugate() ) )
              .divide(f2);
        n22 = n22.multiply(Math.pow(2*Math.PI*freqs[i],2));
        
        Complex n33 = 
            p33.subtract(
                p23.multiply( p13.conjugate() ).divide( p21 ) )
              .divide(f3);
        n33 = n33.multiply(Math.pow(2*Math.PI*freqs[i],2));
        
        // now get magnitude and convert to dB
        double plot1 = 10*Math.log10( n11.abs() );
        double plot2 = 10*Math.log10( n22.abs() );
        double plot3 = 10*Math.log10( n33.abs() );
        if (Math.abs(plot1) != Double.POSITIVE_INFINITY) {
          noiseSeriesArr[0].add(1/freqs[i], plot1);
        }
        if (Math.abs(plot2) != Double.POSITIVE_INFINITY) {
          noiseSeriesArr[1].add(1/freqs[i], plot2);
        }
        if (Math.abs(plot3) != Double.POSITIVE_INFINITY) {
          noiseSeriesArr[2].add(1/freqs[i], plot3);
        }
    }
    
    for (XYSeries noiseSeries : noiseSeriesArr) {
      plottable.addSeries(noiseSeries);
    }
    
    plottable.addSeries( getLowNoiseModel() );
    
    return plottable;
  }
  
  /**
   * Root funtion for calculating crosspower. Gets spectral calculation of data
   * from inputted data series by calling the spectralCalc function, and then
   * applies the provided responses to that result. This is the Power Spectral
   * Density of the inputted data if both sets are the same.
   * @param data1 First data series
   * @param data2 Second data series
   * @param ir1 Response of instrument producing first series
   * @param ir2 Response of instrument producing second series
   * @return Data structure containing the crosspower result of the two data
   * sets as a complex array and the frequencies matched to them in a double
   * array. 
   */
  public FFTStruct crossPower(DataBlock data1, DataBlock data2,
      InstrumentResponse ir1, InstrumentResponse ir2) {
    
    FFTStruct selfPSD = spectralCalc(data1, data2);
    Complex[] results = selfPSD.getFFT();
    double[] freqs = selfPSD.getFreqs();
    Complex[] out = new Complex[freqs.length];
    Complex[] freqRespd1 = ir1.applyResponseToInput(freqs);
    Complex[] freqRespd2 = ir2.applyResponseToInput(freqs);
    freqs = selfPSD.getFreqs();
    
    for (int j = 0; j < freqs.length; ++j) {
      Complex respMagnitude = 
          freqRespd1[j].multiply( freqRespd2[j].conjugate() );
      if (respMagnitude.abs() == 0) {
        respMagnitude = new Complex(Double.MIN_VALUE, 0);
      }
      
      out[j] = results[j].divide(respMagnitude);
      
    }
    
    return new FFTStruct(out, freqs);
    
  }
  
  // TODO: move these signal processing functions into their own
  // class -- or possibly even a new package that can be easily libraried?
  // along with the DataSeriesHelper code
  
  /**
   * Helper function to calculate power spectral density / crosspower.
   * Takes in two time series data and produces the windowed FFT over each.
   * The first is multiplied by the complex conjugate of the second.
   * If the two series are the same, this is the PSD of that series. If they
   * are different, this result is the crosspower.
   * The result is smoothed but does not have the frequency response applied,
   * and so does not give a full result -- this is merely a helper function
   * for the crossPower function.
   * @param data1 DataBlock with relevant time series data
   * @param data2 DataBlock with relevant time series data
   * @return A structure with two arrays: an array of Complex numbers 
   * representing the PSD result, and an array of doubles representing the
   * frequencies of the PSD.
   */
  private FFTStruct spectralCalc(DataBlock data1, DataBlock data2) {
    
    // TODO: try to split off some of the windowed FFT calculations
    // in order to improve quality of code
    
    // this is ugly logic here, but this saves us issues with looping
    // and calculating the same data twice
    boolean sameData = data1.getName().equals(data2.getName());
    
    // divide into windows of 1/4, moving up 1/16 of the data at a time
    
    int range = data1.getData().size()/4;
    int slider = range/4;
    
    // period is 1/sample rate in seconds
    // since the interval data is just that multiplied by a large number
    // let's divide it by that large number to get our period
    
    // shouldn't need to worry about a cast here
    double period = 1.0 / DataBlockHelper.ONE_HZ_INTERVAL;
    period *= data1.getInterval();
    
    int padding = 2;
    while (padding < range) {
      padding *= 2;
    }
    
    int singleSide = padding / 2 + 1;
    double deltaFreq = 1. / (padding * period); // TODO: check this line
    
    Complex[] powSpectDens = new Complex[singleSide];
    double wss = 0;
    
    int segsProcessed = 0;
    int rangeStart = 0;
    int rangeEnd = range;
    
    for (int i = 0; i < powSpectDens.length; ++i) {
      powSpectDens[i] = Complex.ZERO;
    }
    
    while ( rangeEnd <= data1.getData().size() ) {
      
      Complex[] fftResult1 = new Complex[singleSide]; // first half of FFT reslt
      Complex[] fftResult2 = null;
      
      if (!sameData) {
        fftResult2 = new Complex[singleSide];
      }
      
      // give us a new list we can modify to get the data of
      List<Number> data1Range = 
          new ArrayList<Number>(
              data1.getData().subList(rangeStart, rangeEnd) );
      List<Number> data2Range = null;
      
      if (!sameData) {
        data2Range = 
            new ArrayList<Number>(
                data2.getData().subList(rangeStart, rangeEnd) );
      }
       
      // double arrays initialized with zeros, set as a power of two for FFT
      // (i.e., effectively pre-padded on initialization)
      double[] toFFT1 = new double[padding];
      double[] toFFT2 = null;
      
      // demean and detrend work in-place on the list
      demean(data1Range);
      detrend(data1Range);
      wss = cosineTaper(data1Range, TAPER_WIDTH);
      // TODO: check that we only need the last value of wss
      
      if (!sameData) {
        demean(data2Range);
        detrend(data2Range);
        wss = cosineTaper(data2Range, TAPER_WIDTH);
        toFFT2 = new double[padding];
      }
      

      for (int i = 0; i < data1Range.size(); ++i) {
        // no point in using arraycopy -- must make sure each Number's a double
        toFFT1[i] = data1Range.get(i).doubleValue();
        if (!sameData) {
          toFFT2[i] = data2Range.get(i).doubleValue();
        }
      }
      
      FastFourierTransformer fft = 
          new FastFourierTransformer(DftNormalization.STANDARD);


      Complex[] frqDomn1 = fft.transform(toFFT1, TransformType.FORWARD);
      // use arraycopy now (as it's fast) to get the first half of the fft
      System.arraycopy(frqDomn1, 0, fftResult1, 0, fftResult1.length);
      
      Complex[] frqDomn2 = null;
      if (toFFT2 != null) {
        frqDomn2 = fft.transform(toFFT2, TransformType.FORWARD);
        System.arraycopy(frqDomn2, 0, fftResult2, 0, fftResult2.length);
      }
      

      
      for (int i = 0; i < singleSide; ++i) {
        
        Complex val1 = fftResult1[i];
        Complex val2 = val1;
        if (fftResult2 != null) {
          val2 = fftResult2[i];
        }
        
        powSpectDens[i] = 
            powSpectDens[i].add( 
                val1.multiply( 
                    val2.conjugate() ) );
      }
      
      ++segsProcessed;
      rangeStart  += slider;
      rangeEnd    += slider;
      
    }
    
    // normalization time!
    
    double psdNormalization = 2.0 * period / padding; // TODO: check this line
    double windowCorrection = wss / (double) range;
    // it only uses the last value of wss, but that was how the original
    // code was
    
    psdNormalization /= windowCorrection;
    psdNormalization /= segsProcessed; // NOTE: divisor here should be 13
    
    double[] frequencies = new double[singleSide];
    
    for (int i = 0; i < singleSide; ++i) {
      powSpectDens[i] = powSpectDens[i].multiply(psdNormalization);
      frequencies[i] = i * deltaFreq;
    }
    
    // do smoothing over neighboring frequencies; values taken from 
    // asl.timeseries' PSD function
    int nSmooth = 11, nHalf = 5;
    Complex[] psdCFSmooth = new Complex[singleSide];
    
    int iw = 0;
    
    for (iw = 0; iw < nHalf; ++iw) {
      psdCFSmooth[iw] = powSpectDens[iw];
    }
    
    // iw should be icenter of nsmooth point window
    for (; iw < singleSide - nHalf; ++iw){
      int k1 = iw - nHalf;
      int k2 = iw + nHalf;
      
      Complex sumC = Complex.ZERO;
      for (int k = k1; k < k2; ++k) {
        sumC = sumC.add(powSpectDens[k]);
      }
      psdCFSmooth[iw] = sumC.divide((double) nSmooth);
    }
    
    // copy remaining into smoothed array
    for (; iw < singleSide; ++iw) {
      psdCFSmooth[iw] = powSpectDens[iw];
    }
    
    return new FFTStruct(psdCFSmooth, frequencies);
    
  }
  
  /**
   * In-place subtraction of mean from each point in an incoming data set.
   * This is a necessary step in calculating the power-spectral density.
   * @param dataSet The data to have the mean removed from.
   */
  public static void demean(List<Number> dataSet) {
    
    // I'm always getting the demeaning tasks, huh?
    
    if(dataSet.size() == 0) return; // shouldn't happen but just in case
    
    double mean = 0.0;
    
    for(Number data : dataSet) {
      mean += data.doubleValue();
    }
    
    mean /= dataSet.size();
    
    for(int i = 0; i < dataSet.size(); ++i) {
      // iterate over index rather than for-each cuz we must replace data
      // ugly syntax because numeric data types are immutable
      dataSet.set(i, dataSet.get(i).doubleValue() - mean);
    }
    
    // test shows this works as in-place method
  }
  
  /**
   * In-place subtraction of trend from each point in an incoming data set.
   * This is a necessary step in calculating the power-spectral density.
   * @param dataSet The data to have the trend removed from.
   */
  public static void detrend(List<Number> dataSet) {
    
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXSqd = 0.0;
    double sumXY = 0.0;
    
    for (int i = 0; i < dataSet.size(); ++i) {
      sumX += (double) i;
      sumXSqd += (double) i * (double) i;
      double value = dataSet.get(i).doubleValue();
      sumXY += value * (double) i;
      sumY += value;
    }
    
    // brackets here so you don't get confused thinking this should be
    // algebraic division (in which case we'd just factor out the size term)
    // 
    
    double del = sumXSqd - ( sumX * sumX / dataSet.size() );
    
    double slope = sumXY - ( sumX * sumY / dataSet.size() );
    slope /= del;
    
    double yOffset = (sumXSqd * sumY) - (sumX * sumXY);
    yOffset /= del * dataSet.size();
    
    for (int i = 0; i < dataSet.size(); ++i) {
      dataSet.set(i, dataSet.get(i).doubleValue() - ( (slope * i) + yOffset) );
    }
    
  }

  /**
   * Calculates and performs an in-place cosine taper on an incoming data set.
   * Used for windowing for performing FFT.
   * @param dataSet The dataset to have the taper applied to.
   * @return Value corresponding to power loss from application of taper.
   */
  public static double cosineTaper(List<Number> dataSet, double taperW) {
    
    double ramp = taperW * dataSet.size();
    double taper;
    double Wss = 0.0; // represents power loss
    
    for (int i = 0; i < ramp; i++) {
      taper = 0.5 * (1.0 - Math.cos( (double) i * Math.PI / ramp) );
      dataSet.set(i, dataSet.get(i).doubleValue() * taper);
      int idx = dataSet.size()-i-1;
      dataSet.set(idx, dataSet.get(idx).doubleValue() * taper );
      Wss += 2.0 * taper * taper;
    }
    
    Wss += ( dataSet.size() - (2*ramp) );
    
    return Wss;
  }
  
  /**
   * Collects the data points in the Peterson new low noise model to be plotted
   * Assumes that there is a text file in the data folder that contains the
   * NLNM data points for given input frequencies.
   * @return Plottable data series representing the NLNM
   */
  public static XYSeries getLowNoiseModel() {
    // TODO: define NLNM as an array or something in a different class
    XYSeries xys = new XYSeries("NLNM");
    try {
      BufferedReader fr = new BufferedReader(
                            new FileReader(
                              new File("data/NLNM.txt") ) );
      String str = fr.readLine();
      while (str != null) {
        String[] values = str.split("\\s+");
        double x = Double.parseDouble(values[0]);
        if (x > 1.0E3) {
          break;
        }
        double y = Double.parseDouble(values[3]);
        
        xys.add(x,y);
        
        str = fr.readLine();
      }
      fr.close();
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return xys;
  }
  
  /**
   * Holds the data returned from a power spectral density calculation
   * (The PSD data (without response correction) and frequencies of the FFT)
   * @author akearns
   *
   */
  public class FFTStruct {
    
    Complex[] transform;
    double[] freqs;
    
    public FFTStruct(Complex[] inPSD, double[] inFreq) {
      transform = inPSD;
      freqs = inFreq;
    }
    
    public double[] getFreqs() {
      return freqs;
    }
    
    public Complex[] getFFT() {
      return transform;
    }
    
  }
  
  
}
