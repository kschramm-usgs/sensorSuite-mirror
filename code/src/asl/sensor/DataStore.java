package asl.sensor;

import java.io.IOException;

import org.jfree.data.xy.XYSeries;

/**
 * Holds the inputted data from miniSEED files both as a simple struct
 * (see DataBlock) and as a plottable format for the DataPanel.
 * Also serves as container for loaded-in instrument response files
 * @author akearns
 *
 */
public class DataStore {

  /**
   * Defines the maximum number of plots to be shown
   */
  final static int FILE_COUNT = 3;
  DataBlock[] dataBlockArray;
  InstrumentResponse[] responses;
  XYSeries[] outToPlots;
  
  long startTime = 0L;
  long endTime = Long.MAX_VALUE;
  
  boolean[] thisBlockIsSet;
  boolean[] thisResponseIsSet;
  
  /**
   * Instantiate the collections, including empty datasets to be sent to
   * charts for plotting (see DataPanel)
   */
  public DataStore() {
   dataBlockArray = new DataBlock[FILE_COUNT];
   responses = new InstrumentResponse[FILE_COUNT];
   outToPlots = new XYSeries[FILE_COUNT];
   thisBlockIsSet = new boolean[FILE_COUNT];
   thisResponseIsSet = new boolean[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     outToPlots[i] = new XYSeries("(EMPTY) " + i);
     thisBlockIsSet[i] = false;
     thisResponseIsSet[i] = false;
   }
  }
  
  public DataStore(DataStore ds, long start, long end) {
    dataBlockArray = new DataBlock[FILE_COUNT];
    responses = new InstrumentResponse[FILE_COUNT];
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( ds.getBlock(i) != null ) {
        // check in case we're only working with two instead of 3 data sets
        dataBlockArray[i] = new DataBlock( ds.getBlock(i) );
        responses[i] = ds.getResponses()[i];
      }
    }
    // TODO: deep copy?
    thisResponseIsSet = ds.responsesAreSet();
    thisBlockIsSet = ds.dataIsSet();
    this.trimAll(start, end);
  }
  
  public boolean[] dataIsSet() {
    return thisBlockIsSet;
  }
  
  /**
   * Checks if first two blocks and matching responses are set
   * @return True if both of the first two sets have loaded miniSEED and RESP
   */
  public boolean firstTwoSet() {
    // get results for indices 0 and 1 only
    for (int i = 0; i < 2; ++i) {
      if (!thisBlockIsSet[i] || !thisResponseIsSet[i]) {
        return false;
      }
    }
    return true;
  }
  
  /**
   * Return a single data block according to the passed index 
   * @param idx Index of datablock, corresponding to data panel plot index
   * @return Timeseries data for corresponing plot
   */
  public DataBlock getBlock(int idx) {
    if (idx < FILE_COUNT) {
      return dataBlockArray[idx];
    }
    throw new IndexOutOfBoundsException();
  }
  
  /**
   * Returns the set of structures used to hold the loaded miniSeed data sets
   * @return An array of DataBlocks (time series and metadata)
   */
  public DataBlock[] getData() {
    return dataBlockArray;
  }
  
  /**
   * Returns the plottable format of the data held in the arrays at 
   * the specified index
   * @param idx Which of this structure's plottable series should be loaded
   * @return The time series data at given index, to be sent to a chart
   */
  public XYSeries getPlotSeries(int idx) {
    return outToPlots[idx];
  }
  
  /**
   * Returns the instrument responses as an array, ordered such that the
   * response at index i is the response associated with the DataBlock at i
   * in the array of DataBlocks
   * @return
   */
  public InstrumentResponse[] getResponses() {
    return responses;
  }
  
  public boolean isPlottable() {
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (!thisBlockIsSet[i] || !thisResponseIsSet[i]) {
        return false;
      }
    }
    
    return true;
}
  
  public boolean[] responsesAreSet() {
    return thisResponseIsSet;
  }
  
  /**
   * Takes a loaded miniSEED data series and converts it to a plottable format
   * @param idx The plot (range 0 to FILE_COUNT) to be given new data
   * @param filepath Full address of file to be loaded in
   * @return The miniSEED data, in a plottable format
   */
  public XYSeries setData(int idx, String filepath) {
    
    SeedLoaderRunner slr = new SeedLoaderRunner(filepath);
    
    Thread t = new Thread(slr);
    t.start();
    try {
      t.join();
    } catch (InterruptedException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    DataBlock xy = slr.getData(); // executes in a new thread?!
    
    dataBlockArray[idx] = xy;
    
    outToPlots[idx] = xy.toXYSeries();
    
    thisBlockIsSet[idx] = true;
    
    return outToPlots[idx];
  }
  
  /**
   * Sets the response of a sensor's dataseries matched by index
   * @param idx Index of plot for which response file matches
   * @param filepath Full address of file to be loaded in
   */
  public void setResponse(int idx, String filepath) {
    try {
      responses[idx] = new InstrumentResponse(filepath);
      thisResponseIsSet[idx] = true;
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }

  /**
   * Trims all data blocks to be within a certain time range.
   * Used for getting a sub-range specified by sliding-bar window.
   * @param start Start time, relative to epoch (nanoseconds)
   * @param end End time, relative to epoch (nanoseconds)
   */
  public void trimAll(long start, long end) {
    for (DataBlock data : dataBlockArray) {
      if (data != null) {
        data.trim(start, end);
      }
    }
  }

  /**
   * Trims this object's data blocks to hold only points in their common range
   * WARNING: assumes each plot has its data taken at the same point in time
   */
  public void trimToCommonTime() {
    // trims the data to only plot the overlapping time of each input set
    
    long lastStartTime = Long.MIN_VALUE;
    long firstEndTime = Long.MAX_VALUE;
    
    // first pass to get the limits of the time data
    for (DataBlock data : dataBlockArray) {
      if (data == null) {
        continue;
      }
      long start = data.getStartTime();
      if (start > lastStartTime) {
        lastStartTime = start;
      }
      long end = start + data.getInterval() * data.size();
      if (end < firstEndTime) {
        firstEndTime = end;
      }
    }
    
    // second pass to trim the data to the limits given
    for (DataBlock data : dataBlockArray) {
      if (data == null) {
        continue;
      }
      data.trim(lastStartTime, firstEndTime);
    }
    
  }
}
