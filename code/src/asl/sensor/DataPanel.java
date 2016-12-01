package asl.sensor;

import java.awt.Dimension;
import java.awt.Insets;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;

public class DataPanel extends JPanel {

  /**
   * 
   */
  private static final long serialVersionUID = -7302813951637543526L;

  final static int FILE_COUNT = 3;
  
  private boolean[] dataSeries = new boolean[FILE_COUNT];
              // TODO: replace with array of actual data
  private JTextArea[] plotters = new JTextArea[FILE_COUNT]; 
              // replace with plotters for data
  
  public DataPanel() {
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    for (int i = 0; i < FILE_COUNT; ++i) {
      JTextArea statusBox = new JTextArea(5,5);
      statusBox.setMargin( new Insets(5,5,5,5) );
      statusBox.setEditable(false);
      JScrollPane statusScrollPane = new JScrollPane(statusBox);
      statusScrollPane.setBorder( new EmptyBorder(5,5,5,5) );

      
      statusBox.setText("Status Box to eventually hold time-series plots\n");
      statusBox.append("This Box will hold plot for file "+i+'\n');
      
      plotters[i] = statusBox;
      this.add(statusScrollPane);
      this.add( Box.createRigidArea( new Dimension(0,5) ) );
    }
    
  }
  
  public void setData(int idx, boolean data) { // TODO: change param to data fmt
    dataSeries[idx] = true;
    plotters[idx].append("File "+idx+" has been set!\n");
  }
  
  public boolean getData(int idx) { // TODO: change return type to data format
    return dataSeries[idx];
  }
  
  
  
}
