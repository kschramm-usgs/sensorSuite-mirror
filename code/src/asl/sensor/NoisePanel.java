package asl.sensor;

import java.awt.Checkbox;
import java.awt.event.ActionEvent;

import javax.swing.BoxLayout;
import javax.swing.JPanel;

public class NoisePanel extends ExperimentPanel {

  private Checkbox freqSpaceBox;
  
  public NoisePanel(ExperimentEnum exp) {
    // create chart, chartPanel, save button & file chooser, 
    super(exp);
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    // chart = populateChart(expType, expResult);
    
    freqSpaceBox = new Checkbox("Use Hz units (requires regen)");
    freqSpaceBox.setState(false);
    
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
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // only actionlistener here
  }
  
  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    
    if (ds.numberFullySet() < 3) {
      displayErrorMessage("INSUFFICIENT DATA LOADED");
      return;
    }
    
    boolean freqSpace = freqSpaceBox.getState();
    
    updateDriver(ds, psd, freqSpace);
    // setting the new chart is enough to update the plots
    
    
  }
  
  private void updateDriver(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    expResult.setData(ds, psd, freqSpace);
    
    chart = populateChart(expResult.getData(), freqSpace);
    
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);
  }
  

  /**
   * Auto-generated serialize ID
   */
  private static final long serialVersionUID = 9018553361096758354L;
  
  

}
