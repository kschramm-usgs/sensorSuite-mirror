package asl.sensor;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;
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
   * Instantiates a noise experiment -- axis titles and (TODO) scales
   */
  public NoiseExperiment() {
    super();
    xAxisTitle = "Period (s)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
  }

  /**
   * Generates power spectral density of each inputted file, and calculates
   * self-noise based on that result
   * The formula for self-noise calculation is (TODO description)
   */
  @Override
  XYSeriesCollection backend(DataStore ds) {
    
    DataBlock[] dataIn = ds.getData();
    
    XYSeriesCollection plottable = new XYSeriesCollection();
    for (DataBlock data : dataIn) {
      if( data == null ||  data.getData().size() == 0 ) {
        return new XYSeriesCollection();
        // don't actually do any plotting until we have data for everything
      }
      
      
      // first, get crosspower
      
      PSDStruct powSpectResult = powerSpectralDensity(data);
      Complex[] density = powSpectResult.getPSD();
      double[] freqs = powSpectResult.getFreqs();
      
      // next we need to get the instrument response for acceleration
      // and remove that from the PSD
      
    }
    
    return null;
  }
  
  // TODO: move these signal processing functions into their own
  // class -- or possibly even a new package that can be easily libraried?
  // along with the DataSeriesHelper code
  
  /**
   * Helper function to calculate power spectral density
   * @param dataIn DataBlock with relevant time series data
   * @return A structure with two arrays: an array of Complex numbers 
   * representing the PSD result, and an array of doubles representing the
   * frequencies of the PSD.
   */
  public PSDStruct powerSpectralDensity(DataBlock dataIn) {
    
    // divide into windows of 1/4, moving up 1/16 of the data at a time
    
    int range = dataIn.getData().size()/4;
    int slider = range/4;
    
    // period is 1/sample rate in seconds
    // since the interval data is just that multiplied by a large number
    // let's delete it by that large number to get our period
    
    // shouldn't need to worry about a cast here
    double period = 1.0 / DataBlockHelper.ONE_HZ_INTERVAL;
    period *= dataIn.getInterval();
    
    int padding = 2;
    while (padding < range) {
      padding *= 2;
    }
    
    int singleSide = padding / 2 + 1;
    double deltaFreq = 1. / (padding * period);
    
    Complex[] powSpectDens = new Complex[singleSide];
    double wss = 0;
    
    int segsProcessed = 0;
    int rangeStart = 0;
    int rangeEnd = range;
    
    for (int i = 0; i < powSpectDens.length; ++i) {
      powSpectDens[i] = Complex.ZERO;
    }
    
    while ( rangeEnd <= dataIn.getData().size() ) {
      
      Complex[] fftResult = new Complex[singleSide]; // first half of FFT reslt
      
      // give us a new list we can modify to get the data of
      List<Number> dataInRange = 
          new ArrayList<Number>(dataIn.getData().subList(rangeStart, rangeEnd));
      
      
      // demean and detrend work in-place on the list
      demean(dataInRange);
      detrend(dataInRange);
      
      wss = cosineTaper(dataInRange, TAPER_WIDTH);
      
      // double arrays initialized with zeros, set as a power of two for FFT
      // (i.e., effectively pre-padded on initialization)
      double[] toFFT = new double[padding];
      for (int i = 0; i < dataInRange.size(); ++i) {
        // no point in using arraycopy -- must make sure each Number's a double
        toFFT[i] = dataInRange.get(i).doubleValue();
      }
      
      FastFourierTransformer fft = 
          new FastFourierTransformer(DftNormalization.STANDARD);


      Complex[] frqDomn = fft.transform(toFFT, TransformType.FORWARD);
      
      // use arraycopy now (as it's fast) to get the first half of the fft
      System.arraycopy(frqDomn, 0, fftResult, 0, fftResult.length);
      
      for (int i = 0; i < singleSide; ++i) {
        powSpectDens[i] = 
            powSpectDens[i].add( 
                fftResult[i].multiply( 
                    fftResult[i].conjugate() ) );
      }
      
      ++segsProcessed;
      rangeStart  += slider;
      rangeEnd    += slider;
      
    }
    
    // normalization time!
    
    double psdNormalization = 2.0 * period / padding;
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
    
    return new PSDStruct(psdCFSmooth, frequencies);
    
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
   * Holds the data returned from a power spectral density calculation
   * (The PSD data (without response correction) and frequencies of the FFT)
   * @author akearns
   *
   */
  public class PSDStruct {
    
    Complex[] PSD;
    double[] freqs;
    
    public PSDStruct(Complex[] inPSD, double[] inFreq) {
      PSD = inPSD;
      freqs = inFreq;
    }
    
    public double[] getFreqs() {
      return freqs;
    }
    
    public Complex[] getPSD() {
      return PSD;
    }
    
  }
  
  
}
