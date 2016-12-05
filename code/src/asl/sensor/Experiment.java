package asl.sensor;

import org.jfree.data.xy.XYSeriesCollection;

public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: throw exception if dataIn not same size as FILE_COUNT?
  
  protected XYSeriesCollection xySeriesData;
  protected String xAxisTitle, yAxisTitle;
  
  public Experiment(XYSeriesCollection dataIn) {
    xySeriesData = backend(dataIn);
  }
  
  public XYSeriesCollection getData(){
    return xySeriesData;
  }
  
  public void setData(XYSeriesCollection tsc) {
    xySeriesData = backend(tsc);
  }
  
  public String getXTitle() {
    return xAxisTitle;
  }
  
  public String getYTitle() {
    return yAxisTitle;
  }
  
  // java don't allow no static methods 'round here
  abstract XYSeriesCollection backend(XYSeriesCollection tsc);
  
}
