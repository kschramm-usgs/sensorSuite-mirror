package asl.sensor;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.JLabel;

public class StepPanel extends ExperimentPanel {

  JFileChooser fc;
  JButton loadCal;
  JLabel fileNameArea;
  
  
  public StepPanel(ExperimentEnum exp) {
    super(exp);
    
    fc = new JFileChooser();
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    loadCal = new JButton("Load calibration data");
    fileNameArea = new JLabel("NO FILE LOADED");
    
    this.add(chartPanel);
    
    JPanel loadPanel = new JPanel();
    loadPanel.setLayout( new BoxLayout(loadPanel, BoxLayout.X_AXIS) );
    loadPanel.add(loadCal);
    loadPanel.add(fileNameArea);
    this.add(loadPanel);
    
  }

  /**
   * 
   */
  private static final long serialVersionUID = 3693391540945130688L;

  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    // TODO Auto-generated method stub

  }

}
