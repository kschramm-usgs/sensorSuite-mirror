package asl.sensor;

import java.awt.Component;
import java.awt.event.ActionEvent;

import javax.swing.BoxLayout;

public class NoisePanel extends ExperimentPanel {

  public NoisePanel(ExperimentEnum exp) {
    // create chart, chartPanel, save button & file chooser, 
    super(exp);
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    // chart = populateChart(expType, expResult);
    
    // do the layout here
    this.add(chartPanel);
    this.add(save);
    save.setAlignmentX(Component.CENTER_ALIGNMENT);
    

  }
  
  /**
   * Self-noise panel's only interactable object is the save button
   */
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // only actionlistener here
  }

  /**
   * 
   */
  private static final long serialVersionUID = 9018553361096758354L;

  @Override
  public void setDataNames(String[] seedFileNames) {
    // nothing to use the filenames for here
    return;
  }
  
  

}
