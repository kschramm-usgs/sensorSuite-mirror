package asl.sensor;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.Box;
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
import javax.swing.filechooser.FileNameExtensionFilter;


// TODO: include a button for saving plots to PDF somehow
// (double-check on the way that the PDF should be laid out

/**
 * Main window of the sensor test program and the program's launcher
 * Handles GUI for getting user-specified files and showing data plots
 * @author akearns
 *
 */
public class MainWindow extends JPanel implements ActionListener {
  
  /**
   * 
   */
  private static final long serialVersionUID = 2866426897343097822L;
  
  private JButton[] seedLoaders  = new JButton[DataStore.FILE_COUNT];
  private JLabel[] seedFileNames = new JLabel[DataStore.FILE_COUNT];
  private JButton[] respLoaders  = new JButton[DataStore.FILE_COUNT];
  private JLabel[] respFileNames = new JLabel[DataStore.FILE_COUNT];
  
  
  private JFileChooser fc; // loads in files based on parameter
  private DataPanel dataBox;
  private JTabbedPane tabbedPane; // holds set of experiment panels
  
  private JButton generate;
  private JButton savePDF;
  
  
  private void resetTabPlots() {
    DataStore ds = dataBox.getData();
    for ( int i = 0; i < tabbedPane.getTabCount(); ++i ) {
      ExperimentPanel ep = (ExperimentPanel) tabbedPane.getComponentAt(i);
      ep.updateData(ds);
      // updating the chartpanel auto-updates display
    }
    savePDF.setEnabled(true);
  }
  
  /**
   * Creates the main window of the program when called
   * (Three main panels: the top panel for displaying the results
   * of sensor tests; the lower panel for displaying plots of raw data from
   * miniSEED files; the side panel for most file-IO operations
   */
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
    
    for (int i = 0; i < seedLoaders.length; i++){
      seedLoaders[i] = new JButton("Load SEED File " + (i+1) );
      seedLoaders[i].addActionListener(this);
      seedFileNames[i] = new JLabel();
      
      respLoaders[i] = new JButton("Load Response " + (i+1) );
      respLoaders[i].addActionListener(this);
      respFileNames[i] = new JLabel();
      
      initFile(seedLoaders[i], seedFileNames[i], buttonPanel);
      initFile(respLoaders[i], respFileNames[i], buttonPanel);
    }

    generate = new JButton("Generate plots");
    generate.setEnabled(false);
    generate.addActionListener(this);
    buttonPanel.add(generate);
    
    savePDF = new JButton("Save display (PNG)");
    savePDF.setEnabled(true); // TODO: change this back?
    savePDF.addActionListener(this);
    buttonPanel.add(savePDF);
    
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
  
  /**
   * Instantiate a button used to load in a file
   * @param button The button that, when clicked, loads a file
   * @param text Filename (defaults to NO FILE LOADED when created)
   * @param panel The (side) panel that holds the button
   */
  private static void initFile(JButton button, JLabel text, JPanel panel){
    button.setAlignmentX(SwingConstants.CENTER);
    
    text.setText("NO FILE LOADED");
    text.setAlignmentX(SwingConstants.CENTER);

    panel.add(button);
    panel.add(text);
  }
  
  /**
   * Starts the program -- instantiate the top-level GUI
   * @param args (Any parameters fed in on command line are currently ignored)
   */
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
  

  /**
   * Loads the main window for the program on launch
   */
  private static void createAndShowGUI() {
    JFrame frame = new JFrame("Sensor Tests");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    
    frame.add( new MainWindow() );
    
    frame.pack();
    frame.setVisible(true);
  }
  
  /**
   * Handles actions when a side-panel button is clicked (file-IO)
   */
  @Override
  public void actionPerformed(ActionEvent e) {
    // TODO: change control flow 
    
    for(int i = 0; i < seedLoaders.length; ++i) {
      JButton seedButton = seedLoaders[i];
      JButton respButton = respLoaders[i];
      if ( e.getSource() == seedButton ) {
        fc.resetChoosableFileFilters();
        int returnVal = fc.showOpenDialog(seedButton);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File file = fc.getSelectedFile();

          dataBox.setData( i, file.getAbsolutePath() );
          seedFileNames[i].setText( file.getName() );
          
        }
      } else if ( e.getSource() == respButton ) {
        fc.resetChoosableFileFilters();
        int returnVal = fc.showOpenDialog(seedButton);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File file = fc.getSelectedFile();
          
          dataBox.setResponse( i, file.getAbsolutePath() );
          respFileNames[i].setText( file.getName() );
        }
      }
      
      if( dataBox.dataIsSet() ) {
        generate.setEnabled(true);
      }
      
    } // end for loop 
    
    if ( e.getSource() == generate ) {
      this.resetTabPlots();
    } else if ( e.getSource() == savePDF ) {
      
      String ext = ".png";
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      int returnVal = fc.showSaveDialog(savePDF);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          plotsToPNG(selFile);
        } catch (IOException e1) {
          // TODO Auto-generated catch block
          e1.printStackTrace();
        }
      }
      
    }
    

  }
  
  /**
   * Handles function to create a PNG image with all currently-displayed plots
   * @param file File (PNG) that image will be saved to
   * @throws IOException
   */
  public void plotsToPNG(File file) throws IOException {
    // using 0s to set image height and width to default values (match window)
    int inHeight = DataPanel.IMAGE_HEIGHT;
    int width = 640;
    BufferedImage inPlot = dataBox.getAsImage(width, inHeight);
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();

    BufferedImage outPlot = ep.getAsImage( width, 480);
    
    // int width = Math.max( inPlot.getWidth(), outPlot.getWidth() );
    // 5px tall buffer used to separate result plot from inputs
    // BufferedImage space = getSpace(width, 0);
    int height = inPlot.getHeight() + outPlot.getHeight();
    
    // System.out.println(space.getHeight());
    
    BufferedImage toFile = new BufferedImage(width, height, 
        BufferedImage.TYPE_INT_ARGB);
    
    Graphics2D combined = toFile.createGraphics();
    combined.drawImage(outPlot, null, 0, 0);
    combined.drawImage( inPlot, null, 0, 
        outPlot.getHeight() );
    combined.dispose();
    
    // for now, it's a png. TODO: write to PDF?
    
    ImageIO.write(toFile,"png",file);
    
  }
  
  public static BufferedImage getSpace(int width, int height) {
    BufferedImage space = new BufferedImage(
        width, 
        height, 
        BufferedImage.TYPE_INT_RGB);
    Graphics2D tmp = space.createGraphics();
    JPanel margin = new JPanel();
    margin.add( Box.createRigidArea( new Dimension(width,height) ) );
    margin.printAll(tmp);
    tmp.dispose();
    
    return space;
  }
  

}
