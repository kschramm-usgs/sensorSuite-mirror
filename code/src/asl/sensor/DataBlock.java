package asl.sensor;

import java.util.ArrayList;
import java.util.List;

import org.jfree.data.xy.XYSeries;

/**
 * Holds the time series and metadata for a miniSEED file loaded in by the user
 * @author akearns
 *
 */
public class DataBlock {
  
  private List<Number> data;
  private long interval;
  private String name;
  private long startTime;
  
  /**
   * Creates a new DataBlock based on the given parameters
   * @param dataIn Time series data, as a list of numeric objects
   * @param intervalIn The interval between samples (1/sample rate in Hz)
   * @param nameIn The name of the data series, taken from file's metadata
   * @param timeIn The start time of the data series
   */
  public DataBlock
        (List<Number> dataIn, long intervalIn, String nameIn, long timeIn) {
    setData(dataIn);
    setInterval(intervalIn);
    name = nameIn;
    setStartTime(timeIn);
  }
  
  /**
   * Creates a copy of a given DataBlock, which has the same parameters
   * @param in The datablock to be copied
   */
  public DataBlock(DataBlock in) {
    setInterval(in.getInterval());
    setData( new ArrayList<Number>( in.getData() ) );
    name = in.getName();
    setStartTime(in.getStartTime());
  }
  
  /**
   * Get the name of the dataset (used as a key in dataseries, etc.)
   * The name is defined by station-location-channel.
   * @return Name of dataset
   */
  public String getName() {
    return name;
  }
  
  /**
   * Trim data to a given range (start time, end time)
   * @param start Start time to trim to in nanoseconds from epoch
   * @param end End time to trim to in nanoseconds from epoch
   */
  public void trim(long start, long end) {
    int startIdx = 0, endIdx = data.size();
    if (startTime < start) {
      long diff = start - startTime;
      // how many data points in the time range we're removing? 
      startIdx = (int) (diff / interval); // (start offset = 0)
    }
    long endTime = startTime + ( interval * data.size() );
    if ( end < endTime ) {
      long diff = endTime - end;
      // diff/interval is number of points from ending index, need to subtract
      // (quick reminder that upper index of sublist method is exclusive)
      endIdx = endIdx - (int) (diff / interval); // (end offset = size)
    }
    if ( startIdx == 0 && endIdx >= data.size() ){
      return;
    }
    data = new ArrayList<Number>( data.subList(startIdx, endIdx) );
    startTime = start;
    
  }
  
  /**
   * Converts this object's time series data into a form plottable by a chart.
   * The format is a pair of data: the time of a sample and that sample's value.
   * @return JFreeChart XYSeries representation of the data
   */
  public XYSeries toXYSeries() {
    
    final int TRUNCATED_INTERVAL = 100; // used to cheaply decimate plots
    // (the data loaded in is unchanged, this just speeds up display)
    
    XYSeries out = new XYSeries(name);
    long thisTime = getStartTime();
    for (int i = 0; i < data.size(); i+=TRUNCATED_INTERVAL) {
      Number point = data.get(i);
      out.add(thisTime, point);
      thisTime += TRUNCATED_INTERVAL*getInterval();
    }
    
    return out;
  }

  /**
   * Gives the start timestamp of the miniSEED data. This is a long compatible
   * with the Java System Library's Date and Calendar objects and expressed
   * as microseconds from the UTC epoch
   * @return When the miniSEED data logging started in microseconds 
   */
  public long getStartTime() {
    return startTime;
  }

  /**
   * Set the start timestamp of the data with a given long, expressed as
   * microseconds from UTC epoch (compatible with Java System Library Date and
   * Calendar objects)
   * @param startTime The start time in microseconds
   */
  public void setStartTime(long startTime) {
    this.startTime = startTime;
  }

  /**
   * Get the interval of the data. The timestamp for a given data point in the
   * block can be calculated by startTime + (index * interval). This
   * @return The time between two samples of data in microseconds
   */
  public long getInterval() {
    return interval;
  }

  /**
   * Used to set the interval of the data (to be used, for example, when the
   * time series has had decimation applied)
   * @param interval The time between two samples of data in microseconds
   */
  public void setInterval(long interval) {
    this.interval = interval;
  }

  /**
   * Returns the time series object as a list of Java numeric types.
   * The underlying data can be any numeric type, likely Doubles or Integers
   * @return A list representing the miniSEED time series
   */
  public List<Number> getData() {
    return data;
  }

  /**
   * Replace the time series data with a new list of Java numeric types.
   * If the data has been resampled somehow, the resample method is preferred,
   * as it requires a new interval to be resampled
   * @param data New time series data, as from a miniSEED file
   */
  public void setData(List<Number> data) {
    this.data = data;
  }
  
  /**
   * Replace the time series data with a new list of Java numeric types.
   * This is the preferred call when data is resampled as it requires a new
   * interval to be specified, so that the old sample rate does not persist.
   * @param data New time series data, such as the result of decimation
   * @param interval The new interval (time between samples in microseconds)
   */
  public void resample(List<Number> data, long interval) {
    this.data = data;
    this.interval = interval;
  }

  /**
   * Alias to function to get size of list in this data block.
   * Equivalent to calling getData().size()
   * @return The size of the list in this structure
   */
  public int size() {
    return data.size();
  }
  
}
