package asl.sensor;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class DataStore {

  final static int FILE_COUNT = 3;
  DataBlock[] dataBlockArray;
  XYSeries[] outToPlots;
  
  public DataStore(){
   dataBlockArray = new DataBlock[FILE_COUNT];
   outToPlots = new XYSeries[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     outToPlots[i] = new XYSeries("(EMPTY) " + i);
   }
  }
  
  public XYSeries setData(int idx, String filepath) {
    DataBlock xy = DataSeriesHelper.getXYSeries(filepath);
    
    dataBlockArray[idx] = xy;
    outToPlots[idx] = xy.toXYSeries();
    System.out.println(outToPlots[idx].getX(0)+","+outToPlots[idx].getY(0));
    return outToPlots[idx];
     
  }
  
  public XYSeriesCollection getData() {
    XYSeriesCollection xysc = new XYSeriesCollection();
    
    for (XYSeries xys : outToPlots) {
      // TimeSeries reduced = TimeSeriesHelper.reduce(ts);
      xysc.addSeries(xys);
    }
    
    return xysc;
  }
  
  public XYSeries getSeries(int idx) {
    return outToPlots[idx];
  }
}
