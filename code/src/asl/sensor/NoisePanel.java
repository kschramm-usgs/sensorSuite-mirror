package asl.sensor;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;

import javax.swing.JCheckBox;
import javax.swing.JPanel;
import javax.swing.SwingWorker;

import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.data.xy.XYSeriesCollection;

public class NoisePanel extends ExperimentPanel {

  /**
   * Auto-generated serialize ID
   */
  private static final long serialVersionUID = 9018553361096758354L;
  
  private JCheckBox freqSpaceBox;
  private NumberAxis freqAxis;
  
  private final int NOISE_PLOT_COUNT = 6;
  // three PSDs, three self-noise calcs
  
  private final Color[] COLORS = {Color.RED, Color.GREEN, Color.BLUE};
  
  private String freqAxisTitle;
  
  public NoisePanel(ExperimentEnum exp) {
    
    // create chart, chartPanel, save button & file chooser, 
    super(exp);
    
    plotTheseInBold = new String[]{"NLNM","NHNM"};
    
    // instantiate local fields
    xAxisTitle = "Period (s)";
    freqAxisTitle = "Frequency (Hz)";
    yAxisTitle = "Power (rel. 1 (m/s^2)^2/Hz)";
    xAxis = new LogarithmicAxis(xAxisTitle);
    freqAxis = new LogarithmicAxis(freqAxisTitle);
    yAxis = new NumberAxis(yAxisTitle);
    yAxis.setAutoRange(true);
    ( (NumberAxis) yAxis).setAutoRangeIncludesZero(false);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    freqAxis.setLabelFont(bold);
    
    freqSpaceBox = new JCheckBox("Use Hz units (requires regen)");
    freqSpaceBox.setSelected(false);
    
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
    this.add(freqSpaceBox, gbc);
    
    gbc.gridx += 1;
    gbc.weightx = 1.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.anchor = GridBagConstraints.CENTER;
    // gbc.gridwidth = GridBagConstraints.REMAINDER;
    this.add(save, gbc);
    
    // add an empty panel as a spacer to keep the save button in the center
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridx += 1;
    gbc.weightx = 0;
    gbc.anchor = GridBagConstraints.WEST;
    JPanel spacer = new JPanel();
    spacer.setPreferredSize( freqSpaceBox.getPreferredSize() );
    this.add(spacer, gbc);
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    // overridden in the event we add more stuff to this panel
    super.actionPerformed(e); // only actionlistener here
  }
  
  @Override
  public ValueAxis getXAxis() {
    
    // true if using Hz units
    if ( freqSpaceBox.isSelected() ) {
        return freqAxis;
    }
    
    return xAxis;
    
  }
  
  @Override
  public String getXTitle() {
    if (freqSpaceBox.isSelected()) {
      return freqAxisTitle;
    }
    return xAxisTitle;
  }

  @Override
  public void updateData(DataStore ds) {
    
    // TODO: replace with try-catch, put this check in the experiment backend?
    if (ds.numberFullySet() < 3) {
      displayErrorMessage("INSUFFICIENT DATA LOADED");
      return;
    }
    
    boolean freqSpace = freqSpaceBox.isSelected();
    
    updateDriver(ds, freqSpace);
    // setting the new chart is enough to update the plots
    
    
  }
  

  private void updateDriver(DataStore ds, boolean freqSpace) {
    
    final DataStore dsImmutable = ds;
    final boolean freqSpaceImmutable = freqSpace;
    
    displayInfoMessage("Calculating data...");
    
    SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {
      @Override
      public Integer doInBackground() {
        expResult.setData(dsImmutable, freqSpaceImmutable);
        
        XYSeriesCollection xysc = expResult.getData();
        
        for (int i = 0; i < NOISE_PLOT_COUNT; ++i) {
          String name = (String) xysc.getSeriesKey(i);
          Color plotColor = COLORS[i%3];
          seriesColorMap.put(name, plotColor);
          seriesDashedMap.put(name, (i >= 3) );
        }
        
        return 0;
      }

      @Override
      public void done() {
        
        displayInfoMessage("Data loaded...drawing chart");
        
        populateChart(expResult.getData());

        chartPanel.setChart(chart);
        chartPanel.setMouseZoomable(true);
      }

    };
    
    worker.execute();
    
  }
  

}
