package asl.sensor;

import java.awt.BasicStroke;
import java.awt.event.ActionEvent;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class GainPanel extends ExperimentPanel implements ChangeListener {

  /**
   * 
   */
  private static final long serialVersionUID = 6697458429989867529L;
  private JSlider leftSlider;
  private JSlider rightSlider;
  private JComboBox<String> firstSeries;
  private JComboBox<String> secondSeries;
  private JButton recalcButton;
  
  private double low, high;
  
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
    secondSeries = new JComboBox<String>();
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      firstSeries.addItem("FILE " + i + " NOT LOADED");
      secondSeries.addItem("FILE " + i + " NOT LOADED");
    }

    firstSeries.setSelectedIndex(0);
    secondSeries.setSelectedIndex(1);
    
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
    this.add(save); save.setAlignmentX(CENTER_ALIGNMENT);
    
  }

  @Override
  public void setDataNames(String[] seedFileNames) {
    // TODO Auto-generated method stub
    firstSeries.removeAllItems();
    secondSeries.removeAllItems();
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      firstSeries.addItem(seedFileNames[i]);
      secondSeries.addItem(seedFileNames[i]);
    }
    firstSeries.setSelectedIndex(0);
    secondSeries.setSelectedIndex(1);
  }

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
      // if either slider is moved, allow recalc gain and sigma values
      recalcButton.setEnabled(true);
      
      // set the positions of the vertical bars
      
      // get the locations of the sliders and the scale factor
      int leftPos = leftSlider.getValue();
      int rightPos = rightSlider.getValue();
      double scale = (high - low)/1000; // slider range is 0 to 1000
      
      // get the values for where the vertical bars should be
      XYPlot xyp = chartPanel.getChart().getXYPlot();
      double lowFreq = Math.pow(10, low + (scale * leftPos) );
      double highFreq = Math.pow(10, low + (scale * rightPos) );
      
      // draw the bars
      xyp.clearDomainMarkers();
      Marker startMarker = new ValueMarker(lowFreq);
      startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      Marker endMarker = new ValueMarker(highFreq);
      endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      xyp.addDomainMarker(startMarker);
      xyp.addDomainMarker(endMarker);
    }
  }
  
  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    
    int idx1 = firstSeries.getSelectedIndex();
    int idx2 = secondSeries.getSelectedIndex();
    
    // have to make sure the plotting indices get properly set before
    // running the backend, so that we don't add more plots than we need
    expResult.setPlottingIndices(new int[] {idx1, idx2});
    
    super.updateData(ds, psd, true);
    // setting the new chart is enough to update the plots
    
    // now to update all the unique data (ratio, sigma, and the vertical bars)
    double[] meanAndSDev = ( (GainExperiment) expResult).getStats();
    double[] highAndLowFreqs = 
        ( (GainExperiment) expResult).getFreqBoundaries();
    
    double mean = meanAndSDev[0];
    double sDev = meanAndSDev[1];
    
    // set the vertical bars for the plot
    // since we're looking at frequencies, convert to periods
    double lowFreq = highAndLowFreqs[0];
    double highFreq = highAndLowFreqs[1];
    
    double lowPrd = 1. / highFreq;
    double highPrd = 1. / lowFreq;
    
    XYPlot xyp = chartPanel.getChart().getXYPlot();
    xyp.clearDomainMarkers();
    Marker startMarker = new ValueMarker(lowPrd);
    startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    Marker endMarker = new ValueMarker(highPrd);
    endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    xyp.addDomainMarker(startMarker);
    xyp.addDomainMarker(endMarker);
    
    // set slider values to match where the vertical bars are 
    XYSeriesCollection xysc = (XYSeriesCollection) xyp.getDataset();
    XYSeries xys = xysc.getSeries(0);
    if ( xysc.getSeriesKey(0).equals("NLNM") ) {
      xys = xysc.getSeries(1);
    }
    
    // since intervals of incoming data set, should be the same for both inputs
    low = Math.log10( xys.getMinX() ); // value when slider is 0
    high = Math.log10( xys.getMaxX() ); // value when slider is 1000
    
    double scale = (high - low) / 1000; // 1000 <- max value of slider
    
    // refer to http://stackoverflow.com/questions/846221/logarithmic-slider
    int leftValue = (int) ( ( Math.log10(lowPrd) - low ) / scale);
    int rightValue = (int) ( ( Math.log10(highPrd) - low ) / scale);
    
    leftSlider.setEnabled(true);
    rightSlider.setEnabled(true);
    
    leftSlider.setValue(leftValue);
    rightSlider.setValue(rightValue);
    
    // now to add the box with the extra data
    TextTitle result = new TextTitle();
    String temp = "ratio: "+ mean + "\t" + "sigma: " + sDev;
    result.setText(temp);
    
    chartPanel.getChart().addSubtitle(result);
    
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // saving?
    
    if ( e.getSource() == recalcButton ) {
      // we need to get the values of the sliders again
      // TODO: separate out value generation to eliminate duplicated code
      
      int leftPos = leftSlider.getValue();
      int rightPos = rightSlider.getValue();
      double scale = (high - low)/1000; // slider range is 0 to 1000
      
      // get the values for where the vertical bars should be
      // since plot is in seconds (period), need 1/period to get frequency
      double lowFreq = 1. / Math.pow(10, low + (scale * leftPos) );
      double highFreq = 1. / Math.pow(10, low + (scale * rightPos) );
      
      double[] meanAndStdDev = 
          ((GainExperiment) expResult).getStatsFromFreqs(lowFreq, highFreq);
      
      double mean = meanAndStdDev[0];
      double sDev = meanAndStdDev[1];
      
      // TODO: separate out title generation to eliminate duplicated code
      TextTitle result = new TextTitle();
      String temp = "ratio: "+ mean + "\t" + "sigma: " + sDev;
      result.setText(temp);
      chartPanel.getChart().clearSubtitles();
      chartPanel.getChart().addSubtitle(result);
      
      recalcButton.setEnabled(false);
    }
    
    // next need to deal with regenerate ratio/sigma button
  }

}
