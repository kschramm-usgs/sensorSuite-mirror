package asl.sensor;

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
  
  /**
   * Instantiate the collections, including empty datasets to be sent to
   * charts for plotting (see DataPanel)
   */
  public DataStore(){
   dataBlockArray = new DataBlock[FILE_COUNT];
   responses = new InstrumentResponse[FILE_COUNT];
   outToPlots = new XYSeries[FILE_COUNT];
   for (int i = 0; i < FILE_COUNT; ++i) {
     outToPlots[i] = new XYSeries("(EMPTY) " + i);
   }
  }
  
  /**
   * Takes a loaded miniSEED data series and converts it to a plottable format
   * @param idx The plot (range 0 to FILE_COUNT) to be given new data
   * @param filepath Full address of file to be loaded in
   * @return The miniSEED data, in a plottable format
   */
  public XYSeries setData(int idx, String filepath) {
    DataBlock xy = DataBlockHelper.getXYSeries(filepath);
    
    dataBlockArray[idx] = xy;
    outToPlots[idx] = xy.toXYSeries();
    System.out.println(outToPlots[idx].getX(0)+","+outToPlots[idx].getY(0));
    return outToPlots[idx];
     
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
}
