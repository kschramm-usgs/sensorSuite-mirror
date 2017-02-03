package asl.sensor;

import java.awt.Checkbox;
import java.awt.Font;
import java.awt.event.ActionEvent;

import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.SwingWorker;

import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;

public class NoisePanel extends ExperimentPanel {

  private Checkbox freqSpaceBox;
  
  private NumberAxis freqAxis;
  private String freqAxisTitle;
  
  public NoisePanel(ExperimentEnum exp) {
    // create chart, chartPanel, save button & file chooser, 
    super(exp);
    
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
    
    freqSpaceBox = new Checkbox("Use Hz units (requires regen)");
    freqSpaceBox.setState(false);
    
    applyAxesToChart();
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    // chart = populateChart(expType, expResult);
    

    
    JPanel optionsPanel = new JPanel();
    optionsPanel.setLayout( new BoxLayout(optionsPanel, BoxLayout.Y_AXIS) );
    optionsPanel.add(freqSpaceBox);
    optionsPanel.add(save);
    save.setAlignmentX(CENTER_ALIGNMENT);

    
    // do the layout here
    this.add(chartPanel);
    this.add(optionsPanel);

    
  }
  
  @Override
  public ValueAxis getXAxis() {
    
    // true if using Hz units
    if ( freqSpaceBox.getState() ) {
        return freqAxis;
    }
    
    return xAxis;
    
  }
  
  @Override
  public String getXTitle() {
    if (freqSpaceBox.getState()) {
      return freqAxisTitle;
    }
    return xAxisTitle;
  }
  
  @Override
  public String[] seriesToDrawBold() {
    return new String[]{"NLNM"};
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    // overridden in the event we add more stuff to this panel
    super.actionPerformed(e); // only actionlistener here
  }
  
  @Override
  public void updateData(DataStore ds) {
    
    // TODO: replace with try-catch, put this check in the experiment backend?
    if (ds.numberFullySet() < 3) {
      displayErrorMessage("INSUFFICIENT DATA LOADED");
      return;
    }
    
    boolean freqSpace = freqSpaceBox.getState();
    
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
        
        return 0;
      }

      @Override
      public void done() {
        
        displayInfoMessage("Data loaded...drawing chart");
        
        populateChart(expResult.getData());

        chartPanel.setChart(chart);
        chartPanel.setMouseZoomable(false);
      }

    };
    
    new Thread(worker).start();
    return;
    
  }
  

  /**
   * Auto-generated serialize ID
   */
  private static final long serialVersionUID = 9018553361096758354L;
  
  

}
