package asl.sensor;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;

import org.jfree.data.xy.XYSeriesCollection;


public class MainWindow extends JPanel implements ActionListener {
  
  /**
   * 
   */
  private static final long serialVersionUID = 2866426897343097822L;
  private JButton[] fileButtons     = new JButton[DataStore.FILE_COUNT];
  private JLabel[] filenameBoxes = new JLabel[DataStore.FILE_COUNT];
  
  private JFileChooser fc; // loads in files based on parameter
  private DataPanel dataBox;
  private JTabbedPane tabbedPane; // holds set of experiment panels
  
  
  private void resetTabPlots() {
    DataBlock[] tsc = dataBox.getData();
    for ( int i = 0; i < tabbedPane.getTabCount(); ++i ) {
      ExperimentPanel ep = (ExperimentPanel) tabbedPane.getComponentAt(i);
      ep.updateData(tsc);
      // updating the chartpanel auto-updates display
    }
  }
  
  public MainWindow() {
    
    super( new BorderLayout() );
    
    dataBox = new DataPanel();
    
    fc = new JFileChooser();
       
    // each pane will correspond to a plot which gets a test from
    // a test factory; this will return the test corresponding to the plot type
    // which is determined based on an enum of tests
    
    tabbedPane = new JTabbedPane();
    
    for( ExperimentEnum exp : ExperimentEnum.values() ){
      JPanel tab = new ExperimentPanel(exp);
      tab.setLayout( new BoxLayout(tab, BoxLayout.Y_AXIS) );
      tabbedPane.addTab( exp.getName(), tab );
    }
    
    tabbedPane.setBorder( new EmptyBorder(5,0,0,0) );
    
    JPanel buttonPanel = new JPanel();
    buttonPanel.setPreferredSize(new Dimension(100,0));
    buttonPanel.setLayout( new BoxLayout(buttonPanel, BoxLayout.Y_AXIS) );
    
    for (int i = 0; i < fileButtons.length; i++){
      fileButtons[i] = new JButton("Load File " + i);
      fileButtons[i].addActionListener(this);
      filenameBoxes[i] = new JLabel();
      
      initFile(fileButtons[i], filenameBoxes[i], buttonPanel);
    }

    buttonPanel.setBorder( new EmptyBorder(5,5,5,5) );
    
    //holds everything except the side panel used for file IO stuff
    JPanel temp = new JPanel();
    temp.setLayout( new BoxLayout(temp, BoxLayout.Y_AXIS) );
    temp.add(tabbedPane);
    // temp.add(save);
    temp.add(dataBox);
    temp.setBorder( new EmptyBorder(5,5,5,5) );

    JSplitPane splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
    splitPane.setLeftComponent(temp);
    splitPane.setRightComponent(buttonPanel);
    splitPane.setResizeWeight(1.0);
    this.add(splitPane);

  }
  
  private static void initFile(JButton button, JLabel text, JPanel panel){
    button.setAlignmentX(SwingConstants.CENTER);
    
    text.setText("NO FILE LOADED");
    text.setAlignmentX(SwingConstants.CENTER);

    panel.add(button);
    panel.add(text);
  }
  
  public static void main(String[] args) {
    //Schedule a job for the event dispatch thread:
    //creating and showing this application's GUI.
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        //Turn off metal's use of bold fonts
        // UIManager.put("swing.boldMetal", Boolean.FALSE); 
        createAndShowGUI();
      }
    });
    
  }
  

  private static void createAndShowGUI() {
    JFrame frame = new JFrame("Sensor Tests");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    
    frame.add( new MainWindow() );
    
    frame.pack();
    frame.setVisible(true);
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    // TODO: add more checks here as we add components
    
    for(int i = 0; i < fileButtons.length; ++i) {
      JButton button = fileButtons[i];
      if ( e.getSource() == button ) {
        
        // TODO: try-catch on set data to prevent premature renaming
        int returnVal = fc.showOpenDialog(button);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File file = fc.getSelectedFile();

          dataBox.setData(i,file.getAbsolutePath());
          filenameBoxes[i].setText(file.getName());
          this.resetTabPlots();
        }
        return;
      }
    }
  }
  

}
