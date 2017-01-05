package asl.sensor;

import java.awt.event.ActionEvent;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

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
  
  public GainPanel(ExperimentEnum exp) {
    // instantiate common components
    super(exp);
    
    // instantiate unique components
    leftSlider = new JSlider(0, 1000, 0);
    leftSlider.addChangeListener(this);
    rightSlider = new JSlider(0, 1000, 1000);
    rightSlider.addChangeListener(this);
    
    recalcButton = new JButton("Recalc over range");
    recalcButton.setEnabled(false);
    recalcButton.addActionListener(this);
    
    firstSeries = new JComboBox<String>();
    secondSeries = new JComboBox<String>();
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      firstSeries.addItem("FILE NOT LOADED");
      secondSeries.addItem("FILE NOT LOADED");
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
    // TODO Auto-generated method stub
    
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // saving?
  }

}
