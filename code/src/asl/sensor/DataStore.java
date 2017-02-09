package asl.sensor;

import java.io.FileNotFoundException;
import java.io.IOException;

import org.jfree.data.xy.XYSeries;

// TODO (BIG, BIG IMPORTANT TODO) implement locking for file loading and
// for PSD calculation
// need to lock for file read-in to prevent race conditions on time boundary
// need to lock on PSD to prevent redundant calculation thereof

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
  private DataBlock[] dataBlockArray;
  private InstrumentResponse[] responses;
  private XYSeries[] outToPlots;
  FFTResult[] powerSpectra; // used to cache PSD results to speed up calculation
  
  private boolean[] thisBlockIsSet;
  private boolean[] thisResponseIsSet;
  private boolean[] thisPSDIsCalculated;
  
  /**
   * Instantiate the collections, including empty datasets to be sent to
   * charts for plotting (see DataPanel)
   */
  public DataStore() {
   dataBlockArray = new DataBlock[FILE_COUNT];
   responses = new InstrumentResponse[FILE_COUNT];
   outToPlots = new XYSeries[FILE_COUNT];
   powerSpectra = new FFTResult[FILE_COUNT];
   thisBlockIsSet = new boolean[FILE_COUNT];
   thisResponseIsSet = new boolean[FILE_COUNT];
   thisPSDIsCalculated = new boolean[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     outToPlots[i] = new XYSeries("(EMPTY) " + i);
     thisBlockIsSet[i] = false;
     thisResponseIsSet[i] = false;
     thisPSDIsCalculated[i] = false;
   }
  }
  
  public DataStore(DataStore ds, long start, long end) {
    
    dataBlockArray = new DataBlock[FILE_COUNT];
    responses = new InstrumentResponse[FILE_COUNT];
    powerSpectra = new FFTResult[FILE_COUNT];
    thisBlockIsSet = new boolean[FILE_COUNT];
    thisResponseIsSet = new boolean[FILE_COUNT];
    thisPSDIsCalculated = new boolean[FILE_COUNT];
    for (int i = 0; i < FILE_COUNT; ++i) {
      boolean[] setBlocks = ds.dataIsSet();
      boolean[] setResps = ds.responsesAreSet();
      if (setBlocks[i]) {
        dataBlockArray[i] = new DataBlock( ds.getBlock(i) );
        thisBlockIsSet[i] = true;
      }
      if (setResps[i]) {
        responses[i] = ds.getResponse(i);
        thisResponseIsSet[i] = true;
      }
      if (setBlocks[i] && setResps[i]) {
        DataBlock db = dataBlockArray[i];
        long blockStart = db.getStartTime();
        long blockEnd = db.getEndTime();
        if (start == blockStart && end == blockEnd) {
          thisPSDIsCalculated[i] = ds.isPSDSet(i);
          if (thisPSDIsCalculated[i]) {
            powerSpectra[i] = ds.getPSD(i);
          }
        }
      } else {
        thisPSDIsCalculated[i] = false;
      }
    }
    this.trimAll(start, end);
  }
  
  public boolean blockIsSet(int idx) {
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
  
  public boolean[] dataIsSet() {
    return thisBlockIsSet;
  }
  
  public FFTResult[] getAllPSDs() {
    for (int i = 0; i < FILE_COUNT; ++i) {
      if (thisBlockIsSet[i] && thisResponseIsSet[i]) {
        getPSD(i); // calc psd if it's not yet done (but we need data in first)
      }
    }
    return powerSpectra;
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
  
  /**
   * Returns the plottable format of the data held in the arrays at 
   * the specified index
   * @param idx Which of this structure's plottable series should be loaded
   * @return The time series data at given index, to be sent to a chart
   */
  public XYSeries getPlotSeries(int idx) {
    return outToPlots[idx];
  }
  
  public FFTResult getPSD(int idx) {
    if (!thisPSDIsCalculated[idx]) {
      // need both a response and source data
      if (!thisBlockIsSet[idx] || !thisResponseIsSet[idx]) {
        throw new RuntimeException("Not enough data loaded in");
      } else {
        // have enough data but need to calculate before returning
        DataBlock db = dataBlockArray[idx];
        InstrumentResponse ir = responses[idx];
        powerSpectra[idx] = FFTResult.crossPower(db, db, ir, ir);
        thisPSDIsCalculated[idx] = true;
      }
    }

    return powerSpectra[idx];
  }
  
  public InstrumentResponse getResponse(int idx) {
    return responses[idx];
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
  
  public boolean isAnythingSet() {
    for (int i = 0; i< FILE_COUNT; ++i) {
      if (thisBlockIsSet[i] || thisResponseIsSet[i]) {
        return true;
      }
    }
    return false;
  }
  
  public boolean isPSDSet(int idx) {
    return thisPSDIsCalculated[idx];
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

  public void removeData(int idx) {
    dataBlockArray[idx] = null;
    responses[idx] = null;
    outToPlots[idx] = null;
    thisBlockIsSet[idx] = false;
    thisPSDIsCalculated[idx] = false;
    thisResponseIsSet[idx] = false;
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
  public void setData(int idx, String filepath, String nameFilter) {
    
    try {
      DataBlock xy = TimeSeriesUtils.getTimeSeries(filepath, nameFilter);
      thisBlockIsSet[idx] = true;
      dataBlockArray[idx] = xy;
      thisPSDIsCalculated[idx] = false;
      powerSpectra[idx] = null;
      outToPlots[idx] = xy.toXYSeries();
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    if (numberOfBlocksLoaded() > 1) {
      // loading in multiple series of data? trim to common time now
      long start = dataBlockArray[idx].getStartTime();
      long end = dataBlockArray[idx].getEndTime();
      
      // there's clearly already another block loaded, let's make sure they
      // actually have an intersecting time range
      for (int i = 0; i < FILE_COUNT; ++i) {
        if (i != idx && thisBlockIsSet[i]) {
          // whole block either comes before or after the data set
          if (end < dataBlockArray[i].getStartTime() || 
              start > dataBlockArray[i].getEndTime() ) {
            thisBlockIsSet[idx] = false;
            outToPlots[idx] = null;
            dataBlockArray[idx] = null;
            throw new RuntimeException("Time range does not intersect");
          }
        }
      }
      
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
      thisPSDIsCalculated[idx] = false;
      powerSpectra[idx] = null;
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  public boolean timeSeriesSet(int idx) {
    return thisBlockIsSet[idx];
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
    for (int i = 0; i < FILE_COUNT; ++i) {
      DataBlock data = dataBlockArray[i];
      if (data == null) {
        continue;
      }
      data.trim(lastStartTime, firstEndTime);
      outToPlots[i] = data.toXYSeries();
    }
    
  }
}
