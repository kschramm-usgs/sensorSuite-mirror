package asl.sensor;

import java.io.*;
import java.util.HashMap;
import java.util.Map;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;


public class MainWindow extends JPanel implements ActionListener {
  
  private JButton loadButton1, loadButton2, loadButton3;
  private JButton saveButton;
  private JTextArea filenameBox1, filenameBox2, filenameBox3;
  private JFileChooser fc;
  private JTextArea statusBox;
  private JTabbedPane tabbedPane = new JTabbedPane();
  
  public MainWindow() {
    
    super( new BorderLayout() );
    
    statusBox = new JTextArea(10,50);
    statusBox.setMargin( new Insets(5,5,5,5) );
    statusBox.setEditable(false);
    JScrollPane statusScrollPane = new JScrollPane(statusBox);
    statusScrollPane.setBorder( new EmptyBorder(5,5,5,5) );

    
    statusBox.setText("Status Box that will eventually hold time-series plots");
    

    
    fc = new JFileChooser();
       
    // each pane will correspond to a plot which gets a test from
    // a test factory; this will return the test corresponding to the plot type
    // which is determined based on an enum of tests
    
    tabbedPane = new JTabbedPane();
    
    for( Experiment exp : Experiment.values() ){
      JPanel tab = new JPanel();
      
      tab.setLayout( new BoxLayout(tab, BoxLayout.Y_AXIS) );
      JTextArea txa = new JTextArea(10,50);
      txa.setMargin( new Insets(5,5,5,5) );
      txa.setEditable(false);
      tab.add(txa);
      
      JButton save = new JButton("Save Plot");
      save.addActionListener(this);
      tab.add(save);
      
      tabbedPane.addTab( exp.getName(), tab );
      txa.append( "GUI PLOT AREA WIP FOR " + exp.getName() );
    }
    
    tabbedPane.setBorder( new EmptyBorder(5,0,0,0) );
    
    JPanel buttonPanel = new JPanel();
    buttonPanel.setPreferredSize(new Dimension(200,0));
    buttonPanel.setLayout( new BoxLayout(buttonPanel, BoxLayout.Y_AXIS) );
    
    loadButton1 = new JButton("Open file 1");
    loadButton1.addActionListener(this);
    filenameBox1 = new JTextArea(1,10);
    
    initFile(loadButton1, filenameBox1, buttonPanel);
    
    loadButton2 = new JButton("Open file 2");
    loadButton2.addActionListener(this);
    filenameBox2 = new JTextArea(1,10);
    
    initFile(loadButton2, filenameBox2, buttonPanel);
    
    loadButton3 = new JButton("Open file 3");
    loadButton3.addActionListener(this);
    filenameBox3 = new JTextArea(1,10);
    
    initFile(loadButton3, filenameBox3, buttonPanel);

    buttonPanel.setBorder( new EmptyBorder(5,5,5,5) );
    
    //holds everything except the side panel used for file IO stuff
    JPanel temp = new JPanel();
    temp.setLayout( new BoxLayout(temp, BoxLayout.Y_AXIS) );
    temp.add(tabbedPane);
    // temp.add(save);
    temp.add(statusScrollPane);
    temp.setBorder( new EmptyBorder(5,5,5,5) );

    add(temp);
    add(buttonPanel, BorderLayout.EAST);

  }
  
  private static void initFile(JButton button, JTextArea text, JPanel panel){
    JScrollPane filenameScroll;
    
    text.setEditable(false);
    text.setMargin( new Insets(5,5,5,5) );
    text.setText("NO FILE LOADED");
    filenameScroll = new JScrollPane(text);

    panel.add(button);
    panel.add(filenameScroll);
  }
  
  public static void main(String[] args) {
    //Schedule a job for the event dispatch thread:
    //creating and showing this application's GUI.
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        //Turn off metal's use of bold fonts
        UIManager.put("swing.boldMetal", Boolean.FALSE); 
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
    
    if ( e.getSource() == loadButton1 ) {
      // originally used "MainWindow.this" as parameter
      int returnVal = fc.showOpenDialog(loadButton1);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         File file = fc.getSelectedFile();
         filenameBox1.setText(file.getName());
      }
    } else if ( e.getSource() == loadButton2 ) {
      // originally used "MainWindow.this" as parameter
      int returnVal = fc.showOpenDialog(MainWindow.this);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         File file = fc.getSelectedFile();
         filenameBox2.setText(file.getName());
      }
    } else if ( e.getSource() == loadButton3 ) {
      // originally used "MainWindow.this" as parameter
      int returnVal = fc.showOpenDialog(loadButton3);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
         File file = fc.getSelectedFile();
         filenameBox3.setText(file.getName());
      }
    }
  }
  

}
