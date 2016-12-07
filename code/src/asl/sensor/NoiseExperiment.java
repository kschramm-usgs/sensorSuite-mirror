package asl.sensor;

import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeriesCollection;

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
  XYSeriesCollection backend(DataBlock[] dataIn) {
    // TODO Auto-generated method stub
    for (DataBlock data : dataIn) {
      if( data == null ||  data.getData().size() == 0 ) {
        return new XYSeriesCollection();
        // don't actually do any plotting until we have data for everything
      }
    }
    
    return null;
  }
  
  public static void calcCrossPower(DataBlock dataIn) {
    
    
    
  }
  
  /**
   * In-place subtraction of mean from each point in an incoming data set.
   * This is a necessary step in calculating the power-spectral density.
   * @param dataSet The data to have the mean removed from.
   */
  public static void demean(List<Number> dataSet) {
    
    // I'm always getting the demeaning tasks, huh?
    
    // TODO: test that the dataset formed by this is still valid
    
    if(dataSet.size() == 0) return; // shouldn't happen but just in case
    
    double mean = 0.0;
    
    for(Number data : dataSet) {
      mean += data.doubleValue();
    }
    
    mean /= dataSet.size();
    
    for(int i = 0; i < dataSet.size(); ++i) {
      // iterate over index rather than for-each cuz we must replace data
      dataSet.set(i, dataSet.get(i).doubleValue() - mean);
    }
    
    // test will show if we need to return or not or if this works in-place
  }
  
  /**
   * In-place subtraction of trend from each point in an incoming data set.
   * This is a necessary step in calculating the power-spectral density.
   * @param dataSet The data to have the trend removed from.
   */
  public static void detrend(List<Number> dataSet) {
    
    // TODO: also test this
    
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXSqd = 0.0;
    double sumXY = 0.0;
    
    for (int i = 0; i < dataSet.size(); ++i) {
      sumX += (double) i;
      sumXSqd += (double) i * i;
      double value = dataSet.get(i).doubleValue();
      sumXY += value * i;
      sumY += value;
    }
    
    double del = sumXSqd - (sumX * sumX);
    
    double slope = (sumXY - (sumX * sumY) ) / (del);
    
    double yOffset = ( (sumXSqd * sumY) - (sumX * sumXY) ) / (del);
    
    for (int i = 0; i < dataSet.size(); ++i) {
      dataSet.set(i, dataSet.get(i).doubleValue() - ( (slope * i) + yOffset) );
    }
    
  }

  /**
   * Calculates and performs an in-place cosine taper on an incoming data set.
   * @param dataSet The dataset to have the taper applied to.
   * @return Value corresponding to power loss from application of taper.
   */
  public static double cosineTaper(List<Number> dataSet) {
    
    double ramp = TAPER_WIDTH * dataSet.size();
    double taper;
    double Wss = 0.0; // represents power loss
    
    for (int i = 0; i < ramp; i++) {
      taper = 0.5 * (1.0 - Math.cos( (double) i * Math.PI / ramp) );
      dataSet.set(i, dataSet.get(i).doubleValue() * taper);
      int idx = dataSet.size()-i-1;
      dataSet.set(idx, dataSet.get(idx).doubleValue() * taper );
      Wss += 2.0 * taper * taper;
    }
    
    return Wss;
  }
  
  
  
  
}
