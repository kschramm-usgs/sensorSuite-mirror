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
      
      setTitle(mean, sDev);
      
      recalcButton.setEnabled(false);
      
      return;
    } 
    if ( e.getSource() == firstSeries ) {
      
      // if we got here from removing the items from the list
      // (which happens when we load in new data)
      // don't do anything
      if (firstSeries.getItemCount() == 0) {
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
      if (secondSeries.getItemCount() == 0) {
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

  public double mapSliderToPeriod(int position) {
    double scale = (high - low)/1000; // slider range is 0 to 1000
    return Math.pow(10, low + (scale * position) );
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
  
  
  private void setTitle(double mean, double sDev) {
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

  private void setDomainMarkers(double lowPrd, double highPrd, XYPlot xyp) {
    xyp.clearDomainMarkers();
    Marker startMarker = new ValueMarker( lowPrd );
    startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    Marker endMarker = new ValueMarker( highPrd );
    endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
    xyp.addDomainMarker(startMarker);
    xyp.addDomainMarker(endMarker);
  }
  
  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    
    expResult.setData(ds, psd, freqSpace);
    
    updateDataDriver();
  }
  
  private void updateDataDriver() {
      
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

    // recall: we're plotting in period space (secs)
    double lowPrd = xys.getMinX();
    double highPrd = xys.getMaxX();
    
    setDomainMarkers(lowPrd, highPrd, xyp);
    
    // since intervals of incoming data set, should be the same for both inputs
    low = Math.log10( lowPrd ); // value when slider is 0
    high = Math.log10( highPrd ); // value when slider is 1000
    
    
    leftSlider.setEnabled(true);
    rightSlider.setEnabled(true);
    
    leftSlider.setValue(0);
    rightSlider.setValue(1000);
    
    // get statistics from frequency (convert from period)
    double[] meanAndStdDev = 
        ((GainExperiment) expResult).getStatsFromFreqs(
            idx0, idx1, 1 / lowPrd, 1 / highPrd);
    
    double mean = meanAndStdDev[0];
    double sDev = meanAndStdDev[1];
    
    setTitle(mean, sDev);
  }

}
