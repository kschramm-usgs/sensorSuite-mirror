package asl.sensor;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: include axis scale definitions constructed along with axis titles
  // (i.e., logarithmic, etc.)
  
  protected XYSeriesCollection xySeriesData;
  protected String xAxisTitle, freqAxisTitle, yAxisTitle;
  protected NumberAxis xAxis, freqAxis, yAxis;
  protected boolean freqSpace;
  
  /**
   * Return the plottable data for this experiment
   * @return Plottable data 
   */
  public XYSeriesCollection getData() {
    return xySeriesData;
  }
  
  /**
   * Get the names of data series whose timeseries should be plotted in bold
   * (used for reference data such as the NLNM plot for PSD)
   * @return A list of names as Strings
   */
  public String[] getBoldSeriesNames() {
    return new String[]{};
  }
  
  /**
   * Driver to do data processing on inputted data (calls a concrete backend
   * method which is different for each type of experiment)
   * @param ds Timeseries data to be processed
   * @param freqSpace True if the x-axis should be frequency units (False if it
   * should be units of time, for the period)
   */
  public void setData(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    int idx;
    
    idx = ds.getFirstLoadedBlock();
    
    DataBlock[] dataIn = ds.getData();
    
    xySeriesData = new XYSeriesCollection();
    
    long interval = dataIn[idx].getInterval();
    // int length = dataIn[0].size();
    for (DataBlock data : dataIn) {
      
      if ( data == null) {
        // we can have null data for the third entry in some cases
        continue;
      }
      
      if ( data.getInterval() != interval ) {
        throw new RuntimeException("Interval mismatch on datasets.");
      }
    }
    
    backend(ds, psd, freqSpace);
  }
  
  /**
   * Get the name for the domain axis
   * @return the name of the axis
   */
  public String getXTitle() {
    return xAxisTitle;
  }
  
  /**
   * Get the name for the range axis
   * @return the name of the axis
   */
  public String getYTitle() {
    return yAxisTitle;
  }

  /**
   * Get the name for the domain axis when using units of frequency
   * @return the name of the axis
   */
  public String getFreqTitle() {
    return freqAxisTitle;
  }
  
  /**
   * Get the domain axis for plotting
   * @return the axis for the corresponding chart
   */
  public ValueAxis getXAxis() {
    return xAxis;
  }
  
  /**
   * Get the range axis for plotting
   * @return the axis for the corresponding chart
   */
  public ValueAxis getYAxis() {
    return yAxis;
  }
  
  /**
   * Get the domain axis for plotting when using units of frequency
   * @return the axis for the corresponding chart
   */
  public ValueAxis getFreqAxis() {
    return freqAxis;
  }
  
  public static void addToPlot(
      DataBlock[] dataIn, 
      FFTResult[] psd,
      boolean freqSpace,
      XYSeriesCollection xysc,
      Experiment exp) {
    
    // TODO: throw exception if dataIn and psd are diff lengths
    
    for (int i = 0; i < dataIn.length; ++i) {
      
      if( dataIn[i] == null ) {
        continue;
      }
      
      XYSeries powerSeries = new XYSeries( "PSD " + dataIn[i].getName() );

      Complex[] resultPSD = psd[i].getFFT();
      double[] freqs = psd[i].getFreqs();

      for (int j = 0; j < freqs.length; ++j) {
        if (1/freqs[j] > 1.0E3) {
          continue;
        }

        // TODO: is this right (seems to be)
        Complex temp = resultPSD[j].multiply(Math.pow(2*Math.PI*freqs[j],4));

        if (freqSpace) {
          powerSeries.add( freqs[j], 10*Math.log10( temp.abs() ) );
        } else {
          powerSeries.add( 1/freqs[j], 10*Math.log10( temp.abs() ) );
        }
      }

      xysc.addSeries(powerSeries);

    }
    
  }
  
  /**
   * (Overwritten by concrete experiments with specific operations)
   * @param ds Object containing the raw timeseries data to process
   * @param freqSpace True if the data should be plotted by frequency (False
   * if it should be plotted by period) for domain axis
   * @return Plottable series (pl.) of data for all results generated
   */
  abstract void backend(
      DataStore ds, FFTResult[] psd, boolean freqSpace);
  
}
