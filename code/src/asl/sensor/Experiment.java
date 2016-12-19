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
  protected String xAxisTitle, yAxisTitle;
  protected NumberAxis xAxis, yAxis;
  
  public Experiment() {
    
  }
  
  public XYSeriesCollection getData() {
    return xySeriesData;
  }
  
  public String[] getBoldSeriesNames() {
    return new String[]{};
  }
  
  public void setData(DataStore ds) {
    xySeriesData = backend(ds);
  }
  
  public String getXTitle() {
    return xAxisTitle;
  }
  
  public String getYTitle() {
    return yAxisTitle;
  }
  
  public ValueAxis getXAxis() {
    return xAxis;
  }
  
  public ValueAxis getYAxis() {
    return yAxis;
  }
  
  // java don't allow no static methods 'round here
  abstract XYSeriesCollection backend(DataStore ds);
  
}
