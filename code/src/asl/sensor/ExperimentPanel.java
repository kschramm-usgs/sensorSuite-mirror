package asl.sensor;

import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JTextArea;

public class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;
  // likely need a JFreeChart object of some kind here

  JButton save;
  JTextArea txa; // replace with plot object
  private Experiment exp_; // used to define experiment of each plot object
  
  public ExperimentPanel(Experiment exp){
    
    exp_ = exp;
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    save = new JButton("Save Plot");
    save.addActionListener(this);
    
    txa = new JTextArea(10,50);
    txa.setMargin( new Insets(5,5,5,5) );
    txa.setEditable(false);
    txa.setText("GUI PLOT AREA WIP "+exp.getName());
    
    this.add(txa);
    this.add(save);
    
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    // TODO Auto-generated method stub
    if( e.getSource() == save ) {
      // TODO: save plot to file
      txa.append("\nClicked the save button.");
    }
    
  }
  
}
