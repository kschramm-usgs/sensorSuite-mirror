package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.util.HashMap;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.RandomizedExperiment;
import asl.sensor.experiment.ResponseExperiment;
import asl.sensor.input.DataStore;

/**
 * Panel to display results from a randomized calibration experiment.
 * This includes plots of response magnitude and argument, selectable from a
 * drop-down combo box on the panel.
 * The inclusion of two selectable plots means that overrides are necessary
 * to produce output of both plots when creating a report of the results,
 * and that the typical means of assigning the visible chart cannot be used.
 * @author akearns
 *
 */
public class RandomizedPanel extends ExperimentPanel {

  public static final String MAGNITUDE = ResponseExperiment.MAGNITUDE;
  public static final String ARGUMENT = ResponseExperiment.ARGUMENT;
  private static final Color[] COLOR_LIST = 
      new Color[]{Color.RED, Color.BLUE, Color.GREEN};
  /**
   * 
   */
  private static final long serialVersionUID = -1791709117080520178L;
  /**
   * Static helper method for getting the formatted inset string directly
   * from a RandomizedExperiment
   * @param rnd RandomizedExperiment with data to be extracted
   * @return String format representation of data from the experiment
   */
  public static String getInsetString(RandomizedExperiment rnd) {
    
    List<Complex> fitP = rnd.getFitPoles();
    List<Complex> initP = rnd.getInitialPoles();
    
    StringBuilder sbInit = new StringBuilder();
    StringBuilder sbFit = new StringBuilder();
    sbInit.append("INITIAL POLES: \n");
    sbFit.append("FIT POLES: \n");
    for (int i = 0; i < fitP.size(); ++i) {
      sbInit.append( initP.get(i) );
      sbInit.append("  ");
      sbFit.append( fitP.get(i) );
      sbFit.append("  ");
      // want to fit two to a line
      ++i;
      if ( i < fitP.size() ) {
        sbInit.append( initP.get(i) );
        sbInit.append("  ");
        sbFit.append( fitP.get(i) );
        sbFit.append("  ");
      }
      sbInit.append("\n");
      sbFit.append("\n");
    }
    // remove last newline character
    sbFit.deleteCharAt( sbFit.length() - 1 );
    sbInit.append(sbFit);
    return sbInit.toString();
  }
  private ValueAxis degreeAxis;
  private String degreeAxisTitle;
  private JComboBox<String> plotSelection;
  
  private JCheckBox lowFreqBox;

  private JFreeChart magChart, argChart;

  public RandomizedPanel(ExperimentEnum exp) {
    super(exp);
    
    channelType[0] = "Calibration input";
    channelType[1] = "Calibration output from sensor (RESP required)";
    
    yAxisTitle = "10 * log10( RESP(f) )";
    xAxisTitle = "Frequency (Hz)";
    degreeAxisTitle = "phi(RESP(f))";
    
    xAxis = new LogarithmicAxis(xAxisTitle);
    
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    
    degreeAxis = new NumberAxis(degreeAxisTitle);
    degreeAxis.setAutoRange(true);
    
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    
    lowFreqBox = new JCheckBox("Low frequency calibration");
    lowFreqBox.setSelected(true);
    
    applyAxesToChart(); // now that we've got axes defined
    
    // set the GUI components
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.gridwidth = 3;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);
    
    // place the other UI elements in a single row below the chart
    gbc.gridwidth = 1;
    gbc.weighty = 0.0; gbc.weightx = 0.0;
    gbc.anchor = GridBagConstraints.WEST;
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1; gbc.gridx = 0;
    this.add(lowFreqBox, gbc);
    
    gbc.gridx += 1;
    gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    // gbc.gridwidth = GridBagConstraints.REMAINDER;
    this.add(save, gbc);
    
    // plot selection combo box
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridx += 1;
    gbc.weightx = 0;
    gbc.anchor = GridBagConstraints.WEST;
    plotSelection = new JComboBox<String>();
    plotSelection.addItem(MAGNITUDE);
    plotSelection.addItem(ARGUMENT);
    this.add(plotSelection, gbc);
    plotSelection.addActionListener(this);
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    
    super.actionPerformed(e);
    
    if ( e.getSource() == plotSelection ) {
      if (!set) {
        return;
      }
      
      int idx = plotSelection.getSelectedIndex();
      if (idx == 0) {
        chartPanel.setChart(magChart);
      } else {
        chartPanel.setChart(argChart);
      }
      
      return;
      
    }
    
  }
  
  @Override
  public JFreeChart[] getCharts() {
    return new JFreeChart[]{magChart, argChart};
  }
  
  /**
   * Used to get the text that will populate the inset box for the plots
   * @return String to place in TextTitle
   */
  @Override
  public String getInsetString() {
    RandomizedExperiment rnd = (RandomizedExperiment) expResult;
    return getInsetString(rnd);
  }
  
  @Override
  public String getMetadataString() {
    RandomizedExperiment rnd = (RandomizedExperiment) expResult;
    StringBuilder sb = new StringBuilder();
    sb.append("LOADED RESPONSE:");
    sb.append('\n');
    sb.append( rnd.getResponseName() );
    return sb.toString();
  }
  
  @Override
  public ValueAxis getYAxis() {
    
    if ( null == plotSelection ) {
      return yAxis;
    }
    
    if ( plotSelection.getSelectedItem().equals(MAGNITUDE) ) {
      return yAxis;
    } else {
      return degreeAxis;
    }
  }
  
  @Override
  public int panelsNeeded() {
    return 2;
  }

  @Override
  public void updateData(DataStore ds) {

    seriesColorMap = new HashMap<String, Color>();
    
    RandomizedExperiment respExp = (RandomizedExperiment) expResult;
    respExp.setLowFreq( lowFreqBox.isSelected() );
    expResult.setData(ds);
    
    set = true;
    
    List<XYSeriesCollection> xysc = expResult.getData();
    
    XYSeriesCollection magSeries = xysc.get(0);
    XYSeriesCollection argSeries = xysc.get(1);
    
    for (int i = 0; i < magSeries.getSeriesCount(); ++i) {
      
      Color toColor = COLOR_LIST[i];
      String magName = (String) magSeries.getSeriesKey(i);
      seriesColorMap.put(magName, toColor);
      
      String argName = (String) argSeries.getSeriesKey(i);
      seriesColorMap.put(argName, toColor);
      
    }
    
    Range argRange = argSeries.getRangeBounds(true);
    
    int idx = plotSelection.getSelectedIndex();
    
    String inset = getInsetString();
    TextTitle result = new TextTitle();
    result.setText( inset );
    result.setBackgroundPaint(Color.white);

    argChart = buildChart(argSeries);
    argChart.getXYPlot().setRangeAxis(degreeAxis);
    argChart.getXYPlot().getRangeAxis().setAutoRange(true);
    
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
            RectangleAnchor.TOP_RIGHT);
    
    XYPlot xyp = argChart.getXYPlot();
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    magChart = buildChart(magSeries);
    magChart.getXYPlot().setRangeAxis(yAxis);
    magChart.getXYPlot().getRangeAxis().setAutoRange(true);
    
    xyt = new XYTitleAnnotation(0.98, 0.02, result,
        RectangleAnchor.BOTTOM_RIGHT);
    
    xyp = magChart.getXYPlot();
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);

    if (idx == 0) {
      chart = magChart;
    } else {
      chart = argChart;
      chart.getXYPlot().getRangeAxis().setRange(argRange); 
    }
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(true);
    
  }

}
