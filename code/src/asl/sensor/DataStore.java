package asl.sensor;

import java.io.IOException;
import java.util.Arrays;

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
    thisBlockIsSet = new boolean[FILE_COUNT];
    thisResponseIsSet = new boolean[FILE_COUNT];
    for (int i = 0; i < FILE_COUNT; ++i) {
      boolean[] setBlocks = ds.dataIsSet();
      boolean[] setResps = ds.responsesAreSet();
      if ( setBlocks[i] ) {
        dataBlockArray[i] = new DataBlock( ds.getBlock(i) );
        thisBlockIsSet[i] = true;
      }
      if ( setResps[i] ) {
        responses[i] = ds.getResponse(i);
        thisResponseIsSet[i] = true;
      }
    }
    this.trimAll(start, end);
  }
  
  public boolean[] dataIsSet() {
    return thisBlockIsSet;
  }
  
  /**
   * Return a single data block according to the passed index 
   * @param idx Index of datablock, corresponding to data panel plot index
   * @return Timeseries data for corresponing plot
   */
  public DataBlock getBlock(int idx) {
    return dataBlockArray[idx];
  }
  
  /**
   * Returns the set of structures used to hold the loaded miniSeed data sets
   * @return An array of DataBlocks (time series and metadata)
   */
  public DataBlock[] getData() {
    return dataBlockArray;
  }
  
  public InstrumentResponse getResponse(int idx) {
    return responses[idx];
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
  
  public int numberFullySet() {
    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( bothComponentsSet(i) ) {
        ++loaded;
      }
    }
    return loaded;
  }
  
  public int numberOfBlocksLoaded() {
    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i]){
        ++loaded;
      }
    }
    return loaded;
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
  public void setData(int idx, String filepath) {
    
    DataBlock xy = TimeSeriesUtils.getTimeSeriesWithErrorChecking(filepath);
    thisBlockIsSet[idx] = true;
    dataBlockArray[idx] = xy;
    
    if (numberOfBlocksLoaded() > 1) {
      // loading in multiple series of data? trim to common time now
      trimToCommonTime();
    }
    
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
    
    if ( numberOfBlocksLoaded() <= 1 ) {
      return;
    }
    
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
  
  /**
   * Used to get the first, second, etc. loaded block, whether or not it has
   * a loaded response file as well.
   * Used to find the panel where a step calibration is loaded
   * @param x x-th set of data to get, starting at 1
   * @return The Xth DataBlock in this object that is not null
   */
  public DataBlock getXthLoadedBlock(int x) {
    if (x < 1) {
      throw new IndexOutOfBoundsException("Parameter must be >= 1");
    }
    
    int count = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i]) {
        ++count;
        if (count == x) {
          return dataBlockArray[i];
        }
      }
    }
    
    String errMsg = "Not enough data loaded in (found " + count + ")";
    throw new IndexOutOfBoundsException(errMsg);
  }
  
  /**
   * Used to get the first, second, etc. data set loaded. Used when operations
   * reading in data don't require all the inputs to be loaded.
   * Requires both SEED and RESP to be loaded for this to be valid.
   * @param x x-th set of loaded data to get, starting at 1
   * @return index of the loaded data
   */
  public int getXthFullyLoadedIndex(int x) {
    if (x < 1) {
      throw new IndexOutOfBoundsException("Parameter must be >= 1");
    }
    
    int loaded = 0;
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( bothComponentsSet(i) ) {
        ++loaded;
        if (loaded == x) {
          return i;
        }
      }
    }
    
    String errMsg = "Not enough data loaded in (found " + loaded + ")";
    throw new IndexOutOfBoundsException(errMsg);
  }
  
  public boolean timeSeriesSet(int idx) {
    return thisBlockIsSet[idx];
  }

  /**
   * Checks if both components at a specific index are set
   * @param idx Index of data to check if set or not
   * @return True if a miniseed and response have been both loaded in
   */
  public boolean bothComponentsSet(int idx) {
    return (thisBlockIsSet[idx] && thisResponseIsSet[idx]);
  }
  
  public boolean blockIsSet(int idx) {
    return thisBlockIsSet[idx];
  }
}
