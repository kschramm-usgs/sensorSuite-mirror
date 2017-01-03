package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;

public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: include axis scale definitions constructed along with axis titles
  // (i.e., logarithmic, etc.)
  
  protected XYSeriesCollection xySeriesData;
  protected String xAxisTitle, freqAxisTitle, yAxisTitle;
  protected NumberAxis xAxis, freqAxis, yAxis;
  
  /**
   * Default constructor
   */
  public Experiment() {
    
  }
  
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
  public void setData(DataStore ds, boolean freqSpace) {
    xySeriesData = backend(ds, freqSpace);
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
  
  /**
   * (Overwritten by concrete experiments with specific operations)
   * @param ds Object containing the raw timeseries data to process
   * @param freqSpace True if the data should be plotted by frequency (False
   * if it should be plotted by period) for domain axis
   * @return Plottable series (pl.) of data for all results generated
   */
  abstract XYSeriesCollection backend(DataStore ds, boolean freqSpace);
  
}
