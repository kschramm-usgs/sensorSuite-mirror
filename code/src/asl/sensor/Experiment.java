package asl.sensor;

import org.jfree.data.time.TimeSeriesCollection;

public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: throw exception if dataIn not same size as FILE_COUNT?
  
  protected TimeSeriesCollection timeSeriesData;
  protected String xAxisTitle, yAxisTitle;
  
  public Experiment(TimeSeriesCollection dataIn) {
    timeSeriesData = backend(dataIn);
  }
  
  public TimeSeriesCollection getData(){
    return timeSeriesData;
  }
  
  public void setData(TimeSeriesCollection dataIn) {
    timeSeriesData = backend(dataIn);
  }
  
  public String getXTitle() {
    return xAxisTitle;
  }
  
  public String getYTitle() {
    return yAxisTitle;
  }
  
  // java don't allow no static methods 'round here
  abstract TimeSeriesCollection backend(TimeSeriesCollection dataIn);
  
}
