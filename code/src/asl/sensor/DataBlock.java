package asl.sensor;

import java.util.List;

import org.jfree.data.xy.XYSeries;

public class DataBlock {
  
  List<Number> data;
  long interval;
  String name;
  long startTime;
  
  public DataBlock
        (List<Number> dataIn, long intervalIn, String nameIn, long timeIn) {
    data = dataIn;
    interval = intervalIn;
    name = nameIn;
    startTime = timeIn;
  }
  
  public DataBlock(DataBlock in) {
    interval = in.interval;
    data = in.data;
    name = in.name;
    startTime = in.startTime;
  }
  
  public XYSeries toXYSeries() {
    XYSeries out = new XYSeries(name);
    long thisTime = startTime;
    for (Number point : data) {
      out.add(thisTime, point);
      thisTime += interval;
    }
    
    return out;
  }
  
  
}
