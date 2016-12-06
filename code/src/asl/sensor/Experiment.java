package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: include axis scale definitions constructed along with axis titles
  // (i.e., logarithmic, etc.)
  
  protected XYSeriesCollection xySeriesData;
  protected String xAxisTitle, yAxisTitle;
  
  public Experiment() {
    
  }
  
  public XYSeriesCollection getData(){
    return xySeriesData;
  }
  
  public void setData(DataBlock[] tsc) {
    xySeriesData = backend(tsc);
  }
  
  public String getXTitle() {
    return xAxisTitle;
  }
  
  public String getYTitle() {
    return yAxisTitle;
  }
  
  // java don't allow no static methods 'round here
  abstract XYSeriesCollection backend(DataBlock[] tsc);
  
}
