package asl.sensor;

import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

public class DataStore {

  final static int FILE_COUNT = 3;
  TimeSeries[] timeSeriesArray;
  
  public DataStore(){
   timeSeriesArray = new TimeSeries[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     timeSeriesArray[i] = new TimeSeries("");
   }
  }
  
  public TimeSeries setData(int idx, String filepath) {
    TimeSeries ts = TimeSeriesHelper.getTimeSeries(filepath);
    timeSeriesArray[idx] = ts;
    
    return TimeSeriesHelper.reduce(ts);
     
  }
  
  public TimeSeriesCollection getData() {
    TimeSeriesCollection tsc = new TimeSeriesCollection();
    
    for (TimeSeries ts : timeSeriesArray) {
      tsc.addSeries(TimeSeriesHelper.reduce(ts));
    }
    
    return tsc;
  }
  
  public TimeSeries getSeries(int idx) {
    return TimeSeriesHelper.reduce(timeSeriesArray[idx]);
  }
}
