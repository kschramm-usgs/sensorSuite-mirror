package asl.sensor;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: include axis scale definitions constructed along with axis titles
  // (i.e., logarithmic, etc.)
  
  public static void addToPlot(
      final DataStore ds,
      final boolean freqSpace,
      XYSeriesCollection xysc) {
    
    DataBlock[] dataIn = ds.getData();
    
    for (int i = 0; i < dataIn.length; ++i) {
      
      if( dataIn[i] == null ) {
        continue;
      }
      
      XYSeries powerSeries = new XYSeries( "PSD " + dataIn[i].getName() );

      Complex[] resultPSD = ds.getPSD(i).getFFT();
      double[] freqs = ds.getPSD(i).getFreqs();

      for (int j = 0; j < freqs.length; ++j) {
        if (1/freqs[j] > 1.0E3) {
          continue;
        }

        // TODO: is this right (seems to be)
        Complex temp = resultPSD[j].multiply(Math.pow(2*Math.PI*freqs[j],4));

        if (freqSpace) {
          powerSeries.add( freqs[j], 10*Math.log10( temp.abs() ) );
        } else {
          powerSeries.add( 1/freqs[j], 10*Math.log10( temp.abs() ) );
        }
      }

      xysc.addSeries(powerSeries);

    }
    
  }
  
  protected XYSeriesCollection xySeriesData;
  
  /**
   * (Overwritten by concrete experiments with specific operations)
   * @param ds Object containing the raw timeseries data to process
   * @param freqSpace True if the data should be plotted by frequency (False
   * if it should be plotted by period) for domain axis
   * @return Plottable series (pl.) of data for all results generated
   */
  protected abstract void backend(
      final DataStore ds, final boolean freqSpace);
  
  /**
   * Return the plottable data for this experiment
   * @return Plottable data 
   */
  public XYSeriesCollection getData() {
    return xySeriesData;
  }
  
  /**
   * Driver to do data processing on inputted data (calls a concrete backend
   * method which is different for each type of experiment)
   * @param ds Timeseries data to be processed
   * @param freqSpace True if the x-axis should be frequency units (False if it
   * should be units of time, for the period)
   */
  public void setData(final DataStore ds, final boolean freqSpace) {
    

    final DataBlock db = ds.getXthLoadedBlock(1);

    long interval = db.getInterval();
    
    final DataBlock[] dataIn = ds.getData();
    
    xySeriesData = new XYSeriesCollection();
    
    // int length = dataIn[0].size();
    for (final DataBlock data : dataIn) {
      
      if ( data == null) {
        // we can have null blocks, but can't get interval from a null block
        continue;
      }
      
      if ( data.getInterval() != interval ) {
        throw new RuntimeException("Interval mismatch on datasets.");
      }
    }
    
    backend(ds, freqSpace);
  }
  
}
