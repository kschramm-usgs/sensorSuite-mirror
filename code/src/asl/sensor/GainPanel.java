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
      String out = "FILE NOT LOADED (" + i + ")";
      firstSeries.addItem(out);
      secondSeries.addItem(out);
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
      
      // if either slider is moved, allow recalc gain and sigma values
      recalcButton.setEnabled(true);
      
      // get plot (where we put the vertical bars)
      XYPlot xyp = chartPanel.getChart().getXYPlot();
      
      // convert slider locations to (log-scale) frequency
      int leftPos = leftSlider.getValue();
      double lowFreq = mapSliderToFrequency(leftPos);
      int rightPos = rightSlider.getValue();
      double highFreq = mapSliderToFrequency(rightPos);
      
      // remove old bars and draw the new ones
      xyp.clearDomainMarkers();
      Marker startMarker = new ValueMarker(lowFreq);
      startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      xyp.addDomainMarker(startMarker);
      Marker endMarker = new ValueMarker(highFreq);
      endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      xyp.addDomainMarker(endMarker);
    }
  }
  
  public double mapSliderToFrequency(int position) {
    double scale = (high - low)/1000; // slider range is 0 to 1000
    return Math.pow(10, low + (scale * position) );
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
    
    setTitle(mean, sDev);
    
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // saving?
    
    if ( e.getSource() == recalcButton ) {
      // we need to get the values of the sliders again, convert to frequency
      int leftPos = leftSlider.getValue();
      double lowFreq = mapSliderToFrequency(leftPos);
      int rightPos = rightSlider.getValue();
      double highFreq = mapSliderToFrequency(rightPos);
      
      double[] meanAndStdDev = 
          ((GainExperiment) expResult).getStatsFromFreqs(lowFreq, highFreq);
      
      double mean = meanAndStdDev[0];
      double sDev = meanAndStdDev[1];
      
      setTitle(mean, sDev);
      
      recalcButton.setEnabled(false);
    }
    
  }
  
  private void setTitle(double mean, double sDev) {
    // TODO: use XYTitleAnnotation?
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    String temp = "ratio: "+ mean + "\n" + "sigma: " + sDev;
    result.setText(temp);
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }

}
