package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.event.ActionEvent;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class GainPanel extends ExperimentPanel 
implements ChangeListener {

  
  /**
   * 
   */
  private static final long serialVersionUID = 6697458429989867529L;
  private JSlider leftSlider;
  private JSlider rightSlider;
  private JComboBox<String> firstSeries;
  private JComboBox<String> secondSeries;
  private JButton recalcButton;
  boolean freqSpace = false;
  
  private double low, high;
  
  /**
   * Instantiate the panel, including sliders and stat calc button
   * @param exp
   */
  public GainPanel(ExperimentEnum exp) {
    // instantiate common components
    super(exp);
    
    // instantiate unique components
    leftSlider = new JSlider(0, 1000, 0);
    leftSlider.addChangeListener(this);
    leftSlider.setEnabled(false);
    rightSlider = new JSlider(0, 1000, 1000);
    rightSlider.addChangeListener(this);
    rightSlider.setEnabled(false);
    
    recalcButton = new JButton("Recalc over range");
    recalcButton.setEnabled(false);
    recalcButton.addActionListener(this);
    
    firstSeries = new JComboBox<String>();
    firstSeries.addActionListener(this);
    secondSeries = new JComboBox<String>();
    secondSeries.addActionListener(this);
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      String out = "FILE NOT LOADED (" + i + ")";
      firstSeries.addItem(out);
      secondSeries.addItem(out);
    }

    // we will force these to be disabled to load in first two data sets
    // in all cases -- so we can only load in two inputs and get that value
    firstSeries.setSelectedIndex(0);
    firstSeries.setEnabled(false);
    secondSeries.setSelectedIndex(1);
    secondSeries.setEnabled(false);
    
    // create layout
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    JPanel sliderPanel = new JPanel();
    sliderPanel.setLayout( new BoxLayout(sliderPanel, BoxLayout.X_AXIS) );
    sliderPanel.add(leftSlider);
    sliderPanel.add(recalcButton);
    sliderPanel.add(rightSlider);

    JPanel comboPanel = new JPanel();
    comboPanel.setLayout( new BoxLayout(comboPanel, BoxLayout.X_AXIS) );
    comboPanel.add(firstSeries);
    comboPanel.add(secondSeries);
    
    this.add(chartPanel);
    this.add(sliderPanel);
    this.add(comboPanel);
    // save button instantiated but not displayed on panel
    this.add(save); save.setAlignmentX(CENTER_ALIGNMENT);
    
  }

  /**
   * Replots when a new data source is chosen (inc. check to prevent same data
   * from being plotted twice -- which causes an error with JFreeChart)
   * and recalculates stats when the recalculate button is hit
   */
  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // saving?
    
    int idx0 = firstSeries.getSelectedIndex();
    int idx1 = secondSeries.getSelectedIndex();
    
    if ( e.getSource() == recalcButton ) {
      
      // we need to get the values of the sliders again, convert to frequency
      int leftPos = leftSlider.getValue();
      double lowPrd = mapSliderToPeriod(leftPos);
      int rightPos = rightSlider.getValue();
      double highPrd = mapSliderToPeriod(rightPos);
      
      double[] meanAndStdDev = 
          ((GainExperiment) expResult).getStatsFromFreqs(
              idx0, idx1, 1/lowPrd, 1/highPrd);
      
      double mean = meanAndStdDev[0];
      double sDev = meanAndStdDev[1];
      double refGain = meanAndStdDev[2];
      double calcGain = meanAndStdDev[3];
      
      setTitle(mean, sDev, refGain, calcGain);
      
      recalcButton.setEnabled(false);
      
      return;
    } 
    if ( e.getSource() == firstSeries ) {
      
      // if we got here from removing the items from the list
      // (which happens when we load in new data)
      // don't do anything
      if (firstSeries.getItemCount() < 2) {
        return;
      }
      
      // don't allow the same value for the two indices (plot behaves badly)
      // assume the user is setting firstSeries selection correctly,
      // and thus make sure that the secondSeries doesn't have a collision
      if (idx0 == idx1) {
        secondSeries.setSelectedIndex( 
            (idx1 + 1) % secondSeries.getItemCount() );
      }
      
    } else if ( e.getSource() == secondSeries ) {      
      
      // same as above, do nothing
      if (secondSeries.getItemCount() < 2) {
        return;
      }
      
      // same as with the above, but assume secondSeries selection intentional
      if (idx0 == idx1) {
        firstSeries.setSelectedIndex(
            (idx0 + 1) % firstSeries.getItemCount() );
      }
      
    }
    // now that we have a guarantee of no collision, update data accordingly    
    if ( e.getSource() == firstSeries || e.getSource() == secondSeries) {
      // if we selected a new series to plot, redraw the chart
      if (expResult.getData() != null) {
        updateDataDriver();

      }
    }
    
  }

  /**
   * Converts the slider position to a logarithmic scale matching x-axis values
   * which is the period given in a rate of seconds
   * @param position value of slider
   * @return x-axis value corresponding to that position
   */
  public double mapSliderToPeriod(int position) {
    double scale = (high - low)/1000; // slider range is 0 to 1000
    return Math.pow(10, low + (scale * position) );
  }
  
  /**
   * Converts x-axis value from log scale to linear, to get slider position
   * @param prd period value marking data window boundary
   * @return value of slider (ranges from 0 to 1000)
   */
  public int mapPeriodToSlider(double prd) {
    double scale = (high - low)/1000; // recall slider range is 0 to 1000
    return (int) ( ( Math.log10(prd) - low ) / scale );
  }
  
  /**
   * Used to populate the combo box with the names of inputted data
   */
  @Override
  public void setDataNames(String[] seedFileNames) {
    firstSeries.removeAllItems();
    secondSeries.removeAllItems();
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      String out = seedFileNames[i] + " (" + i + ")";
      firstSeries.addItem(out);
      secondSeries.addItem(out);
    }
    firstSeries.setSelectedIndex(0);
    secondSeries.setSelectedIndex(1);
  }
  
  
  /**
   * Displays the statistic results when the calculate button is hit
   * @param mean Calculated mean value
   * @param sDev Calculated standard deviation value
   */
  private void 
  setTitle(double mean, double sDev, double refGain, double calcGain) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    StringBuilder sb = new StringBuilder();
    sb.append("ratio: ");
    sb.append(mean);
    sb.append("\n");
    sb.append("sigma: ");
    sb.append(sDev);
    sb.append("\n");
    sb.append("ref. gain: ");
    sb.append(refGain);
    sb.append("\n");
    sb.append("** CALCULATED GAIN: ");
    sb.append(calcGain);
    String temp = sb.toString();
    result.setText(temp);
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
    
  /**
   * Used to get boundaries of chart window specified by this panel's slider
   * and to draw the vertical lines matching the location of those bounds
   */
  @Override
  public void stateChanged(ChangeEvent e) {
    
    // enforce slider boundaries
    if ( e.getSource() == leftSlider ) {
      if ( leftSlider.getValue() > rightSlider.getValue() - 10 ) {
        leftSlider.setValue( rightSlider.getValue() - 10 );
        if ( leftSlider.getValue() < 0 ) {
          leftSlider.setValue(0);
          rightSlider.setValue(10);
        }
      }
    } else if ( e.getSource() == rightSlider ) {
      if ( leftSlider.getValue() + 10 > rightSlider.getValue() ) {
        rightSlider.setValue( leftSlider.getValue() + 10 );
        if ( rightSlider.getValue() > 1000 ) {
          rightSlider.setValue(1000);
          leftSlider.setValue(990);
        }
      }
    }
    
    if (e.getSource() == leftSlider || e.getSource() == rightSlider) {
      
      // now we need to redraw the vertical bars
      
      // new slider window means new results on calculation
      recalcButton.setEnabled(true);
      
      // get plot (where we put the vertical bars)
      XYPlot xyp = chartPanel.getChart().getXYPlot();
      
      // clear out annotations to prevent issues with misleading data
      xyp.clearAnnotations();
      
      // convert slider locations to (log-scale) frequency
      int leftPos = leftSlider.getValue();
      double lowPrd = mapSliderToPeriod(leftPos);
      int rightPos = rightSlider.getValue();
      double highPrd = mapSliderToPeriod(rightPos);
      
      // remove old bars and draw the new ones
      setDomainMarkers(lowPrd, highPrd, xyp);
    }
  }

  /**
   * Draws the lines marking the boundaries of the current window
   * @param lowPrd lower x-axis value (period, in seconds)
   * @param highPrd upper x-axis value (period, in seconds)
   * @param xyp plot displayed in this object's chart
   */
  private void setDomainMarkers(double lowPrd, double highPrd, XYPlot xyp) {
    xyp.clearDomainMarkers();
    Marker startMarker = new ValueMarker( lowPrd );
    startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    Marker endMarker = new ValueMarker( highPrd );
    endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    xyp.addDomainMarker(startMarker);
    xyp.addDomainMarker(endMarker);
  }
  
  /**
   * Sets the data and then calls function to update data based on
   * the selected combo box entries
   */
  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    
    if (ds.amountOfDataLoaded() < 2) {
      displayErrorMessage("INSUFFICIENT DATA LOADED");
      return;
    }
    
    expResult.setData(ds, psd, freqSpace);
    // need to have 2 series for relative gain
    if ( ds.amountOfDataLoaded() > 2 ) {
      firstSeries.setEnabled(true);
      secondSeries.setEnabled(true);
    } else {
      firstSeries.setEnabled(false);
      secondSeries.setEnabled(false);
    }
    
    updateDataDriver();
  }
  
  /**
   * Given input data (including time series collection), get only the relevant
   * ones to display based on combo boxes and then do the statistics on those
   * Called when new data is loaded or when the combo box active entries change
   */
  private void updateDataDriver() {
    
    // these should be 0 and 1 since the series cannot be selected currently
    // TODO: add check if 3 data sets are loaded in to allow choosing?
    int idx0 = firstSeries.getSelectedIndex();
    int idx1 = secondSeries.getSelectedIndex();
    
    // have to make sure the plotting indices get properly set before
    // running the backend, so that we don't add more plots than we need
     
    XYSeriesCollection xysc = expResult.getData();
    
    XYSeriesCollection plotXYSC = new XYSeriesCollection();
      
    plotXYSC.addSeries( xysc.getSeries(idx0) );
    plotXYSC.addSeries( xysc.getSeries(idx1) );
    plotXYSC.addSeries( xysc.getSeries("NLNM") );
    
    chart = populateChart(plotXYSC, freqSpace);
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);
    
    // set vertical bars and enable sliders
    XYPlot xyp = chartPanel.getChart().getXYPlot();
    
    XYSeries xys = xysc.getSeries(0);
    if ( xysc.getSeriesKey(0).equals("NLNM") ) {
      xys = xysc.getSeries(1);
    }

    // want to default to octave centered at highest value of fixed freq
    double[] freqRange = 
        ( (GainExperiment) expResult).getOctaveCenteredAtPeak(idx0);
    
    // get the locations (x-axis values) of frequency range in 
    double lowPrd = Math.min(1/freqRange[0], 1/freqRange[1]);
    double highPrd = Math.max(1/freqRange[0], 1/freqRange[1]);
    
    // since intervals of incoming data match, so too limits of plot
    // this is used in mapping scale of slider to x-axis values
    low = Math.log10( xys.getMinX() ); // value when slider is 0
    high = Math.log10( xys.getMaxX() ); // value when slider is 1000
    
    // set the domain to match the boundaries of the octave centered at peak
    setDomainMarkers(lowPrd, highPrd, xyp);
    
    // and now set the sliders to match where that window is
    leftSlider.setEnabled(true);
    rightSlider.setEnabled(true);
    
    int leftSliderValue = mapPeriodToSlider(lowPrd);
    int rightSliderValue = mapPeriodToSlider(highPrd);
    
    leftSlider.setValue(leftSliderValue);
    rightSlider.setValue(rightSliderValue);
    
    // get statistics from frequency (convert from period)
    double[] gainStatistics = 
        ( (GainExperiment) expResult ).getStatsFromFreqs(
            idx0, idx1, freqRange[0], freqRange[1]);
    
    double mean = gainStatistics[0];
    double sDev = gainStatistics[1];
    double refGain = gainStatistics[2];
    double calcGain = gainStatistics[3];
    
    setTitle(mean, sDev, refGain, calcGain);
  }

}
