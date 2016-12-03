package asl.sensor;

import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

public class DataStore {

  final static int FILE_COUNT = 3;
  TimeSeries[] timeSeriesArray;
  TimeSeries[] timeSeriesReduced;
  
  public DataStore(){
   timeSeriesArray = new TimeSeries[FILE_COUNT];
   timeSeriesReduced = new TimeSeries[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     timeSeriesArray[i] = new TimeSeries("");
     timeSeriesReduced[i] = new TimeSeries("");
   }
  }
  
  public TimeSeries setData(int idx, String filepath) {
    TimeSeries ts = TimeSeriesHelper.getTimeSeries(filepath);
    timeSeriesArray[idx] = ts;
    timeSeriesReduced[idx] = TimeSeriesHelper.reduce(ts);
    return timeSeriesReduced[idx];
     
  }
  
  public TimeSeriesCollection getData() {
    TimeSeriesCollection tsc = new TimeSeriesCollection();
    
    for (TimeSeries ts : timeSeriesReduced) {
      // TimeSeries reduced = TimeSeriesHelper.reduce(ts);
      tsc.addSeries(ts);
    }
    
    return tsc;
  }
  
  public TimeSeries getSeries(int idx) {
    return timeSeriesReduced[idx];
  }
}
