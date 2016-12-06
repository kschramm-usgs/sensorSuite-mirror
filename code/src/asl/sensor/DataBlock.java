package asl.sensor;

import java.util.List;

import org.jfree.data.xy.XYSeries;

public class DataBlock {
  
  private List<Number> data;
  private long interval;
  String name;
  private long startTime;
  
  public DataBlock
        (List<Number> dataIn, long intervalIn, String nameIn, long timeIn) {
    setData(dataIn);
    setInterval(intervalIn);
    name = nameIn;
    setStartTime(timeIn);
  }
  
  public DataBlock(DataBlock in) {
    setInterval(in.getInterval());
    setData(in.getData());
    name = in.name;
    setStartTime(in.getStartTime());
  }
  
  public XYSeries toXYSeries() {
    XYSeries out = new XYSeries(name);
    long thisTime = getStartTime();
    for (Number point : getData()) {
      out.add(thisTime, point);
      thisTime += getInterval();
    }
    
    return out;
  }

  public long getStartTime() {
    return startTime;
  }

  public void setStartTime(long startTime) {
    this.startTime = startTime;
  }

  public long getInterval() {
    return interval;
  }

  public void setInterval(long interval) {
    this.interval = interval;
  }

  public List<Number> getData() {
    return data;
  }

  public void setData(List<Number> data) {
    this.data = data;
  }
  
  
}
