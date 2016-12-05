package asl.sensor;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class DataStore {

  final static int FILE_COUNT = 3;
  XYSeries[] xySeriesArray;
  
  public DataStore(){
   xySeriesArray = new XYSeries[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     xySeriesArray[i] = new XYSeries("(EMPTY)" + i);
   }
  }
  
  public XYSeries setData(int idx, String filepath) {
    XYSeries xy = DataSeriesHelper.getXYSeries(filepath);
    xySeriesArray[idx] = xy;
    return xy;
     
  }
  
  public XYSeriesCollection getData() {
    XYSeriesCollection xysc = new XYSeriesCollection();
    
    for (XYSeries xys : xySeriesArray) {
      // TimeSeries reduced = TimeSeriesHelper.reduce(ts);
      xysc.addSeries(xys);
    }
    
    return xysc;
  }
  
  public XYSeries getSeries(int idx) {
    return xySeriesArray[idx];
  }
}
