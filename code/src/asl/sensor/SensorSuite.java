package asl.sensor;

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.imageio.ImageIO;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.border.EmptyBorder;
import javax.swing.filechooser.FileNameExtensionFilter;

// TODO: move file operations to the datapanel object?

/**
 * Main window of the sensor test program and the program's launcher
 * Handles GUI for getting user-specified files and showing data plots
 * @author akearns
 *
 */
public class SensorSuite extends JPanel implements ActionListener {

  /**
   * 
   */
  private static final long serialVersionUID = 2866426897343097822L;

  ExecutorService executor = Executors.newCachedThreadPool();
  
  private JButton[] seedLoaders  = new JButton[DataStore.FILE_COUNT];
  private JTextField[] seedFileNames = new JTextField[DataStore.FILE_COUNT];
  private JButton[] respLoaders  = new JButton[DataStore.FILE_COUNT];
  private JTextField[] respFileNames = new JTextField[DataStore.FILE_COUNT];

  private JFileChooser fc; // loads in files based on parameter
  private DataPanel dataBox;
  private JTabbedPane tabbedPane; // holds set of experiment panels

  private JButton generate; // run all calculations
  private JButton savePDF;
  private JButton clear;
  
  // used to store current directory locations
  private String seedDirectory = "data";
  private String respDirectory = "responses";
  private String saveDirectory = "";


  private void resetTabPlots() {
    
    DataStore ds = dataBox.getData();
    InstrumentResponse[] irs = ds.getResponses();
    
    FFTResult[] powerSpectra = new FFTResult[DataStore.FILE_COUNT];
    
    dataBox.showRegionForGeneration();
    
    for (int i = 0; i < powerSpectra.length; ++i) {
      
      if ( !ds.bothComponentsSet(i) ) {
        continue;
      }
      
      DataBlock data = ds.getBlock(i);
      InstrumentResponse ir = irs[i];
      
      powerSpectra[i] = FFTResult.crossPower(data, data, ir, ir);
    }
    
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
    ep.updateData(ds, powerSpectra);
    
    /*
    for ( int i = 0; i < tabbedPane.getTabCount(); ++i ) {
    
    boolean selectedIsTwoInput = 
        ( (ExperimentPanel) tabbedPane.getSelectedComponent() ).isTwoInput();
    
      ExperimentPanel ep = (ExperimentPanel) tabbedPane.getComponentAt(i);
      if (!allThree) {
        // triggered if only two inputs exist (not all three)
        if ( !ep.isTwoInput() ) {
          // ignore anything that requires more than two inputs then
          continue;
        } else if ( !selectedIsTwoInput ) {
          // display the first panel that actually has data to update with
          tabbedPane.setSelectedComponent(ep);
          selectedIsTwoInput = true;
        }
      }
    
      // now update the panel with the correct data to plot
      ep.updateData(ds, powerSpectra);
      // updating the chartpanel auto-updates display
    }
    */
    savePDF.setEnabled(true);
  }

  /**
   * Creates the main window of the program when called
   * (Three main panels: the top panel for displaying the results
   * of sensor tests; the lower panel for displaying plots of raw data from
   * miniSEED files; the side panel for most file-IO operations
   */
  public SensorSuite() {

    super();
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints c = new GridBagConstraints();
    c.fill = GridBagConstraints.BOTH;
    c.gridx = 0; c.gridy = 0;
    c.weightx = 1.0; c.weighty = 1.0;
    c.gridwidth = 1;
    c.anchor = GridBagConstraints.CENTER;
    dataBox = new DataPanel();

    fc = new JFileChooser();
    
    // each pane will correspond to a plot which gets a test from
    // a test factory; this will return the test corresponding to the plot type
    // which is determined based on an enum of tests

    tabbedPane = new JTabbedPane();

    for( ExperimentEnum exp : ExperimentEnum.values() ){
      JPanel tab = ExperimentPanelFactory.createPanel(exp);
      tab.setLayout( new BoxLayout(tab, BoxLayout.Y_AXIS) );
      tabbedPane.addTab( exp.getName(), tab );
    }
    
    this.add(tabbedPane, c);
    
    // reset for other components
    c.weightx = 1.0;
    c.weighty = 1.0;
    c.fill = GridBagConstraints.BOTH;

    // panel for UI elements to load data, generate full plot, etc.
    JPanel rightPanel = new JPanel();
    
    rightPanel.setLayout( new GridBagLayout() );
    GridBagConstraints rpc = new GridBagConstraints();
    rpc.weightx = 0.0;
    rpc.gridwidth = 1;
    rpc.gridy = 0;
    
    rpc.weighty = 1.0;
    rpc.fill = GridBagConstraints.BOTH;
    rightPanel.add( Box.createVerticalGlue(), rpc);
    
    rpc.gridy += 1;
    
    rpc.weighty = 0.0;
    rpc.fill = GridBagConstraints.NONE;

    for (int i = 0; i < seedLoaders.length; i++){
      seedLoaders[i] = new JButton("Load SEED File " + (i+1) );
      seedLoaders[i].addActionListener(this);
      seedFileNames[i] = new JTextField();

      respLoaders[i] = new JButton("Load Response " + (i+1) );
      respLoaders[i].addActionListener(this);
      respFileNames[i] = new JTextField();

      seedFileNames[i].setEditable(false);
      respFileNames[i].setEditable(false);

      initFile(seedLoaders[i], seedFileNames[i], rightPanel, rpc);
      initFile(respLoaders[i], respFileNames[i], rightPanel, rpc);

      rpc.gridy += 1;
      
      rpc.weighty = 1.0;
      rpc.fill = GridBagConstraints.VERTICAL;
      rightPanel.add( Box.createVerticalGlue(), rpc);
      rpc.gridy += 1;
      
      rpc.weighty = 0.0;
      rpc.fill = GridBagConstraints.NONE;
    }
    

    rpc.weighty = 1.0;
    rpc.fill = GridBagConstraints.VERTICAL;
    rightPanel.add( Box.createVerticalGlue(), rpc );
    rpc.gridy += 1;
    
    rpc.weighty = 0.0;
    rpc.fill = GridBagConstraints.HORIZONTAL;
    rpc.ipady = 10;
    
    // TODO: replace duplicated effects factory-style methods?
    
    generate = new JButton("Generate plots");
    generate.setEnabled(true);
    generate.addActionListener(this);
    rightPanel.add(generate, rpc);
    rpc.gridy += 1;

    savePDF = new JButton("Save display (PNG)");
    savePDF.setEnabled(true); // TODO: change this back?
    savePDF.addActionListener(this);
    rightPanel.add(savePDF, rpc);
    rpc.gridy += 1;
    
    clear = new JButton("Clear data");
    clear.setEnabled(true);
    clear.addActionListener(this);
    rightPanel.add(clear, rpc);
    // TODO: clear out the old data in the result plot too 
    
    // add space between the plots and the file-operation panel
    dataBox.setBorder( new EmptyBorder(0, 0, 0, 5) );
    
    JPanel data = new JPanel();
    data.setLayout( new GridBagLayout() );
    
    GridBagConstraints dgbc = new GridBagConstraints();
    dgbc.weightx = 1.0; dgbc.weighty = 1.0;
    dgbc.gridheight = 1;
    dgbc.gridy = 0; dgbc.gridx = 0;
    dgbc.fill = GridBagConstraints.BOTH;
    
    data.add(dataBox, dgbc);
    
    dgbc.weightx = 0.0;
    dgbc.gridx += 1;
    
    data.add(rightPanel, dgbc);
    data.setBorder( new EmptyBorder(5, 5, 5, 5) );
    
    c.gridy += 1;
    this.add(data, c);
    data.setAlignmentX(SwingConstants.CENTER);

  }

  /**
   * Instantiate a button used to load in a file
   * @param button The button that, when clicked, loads a file
   * @param text Filename (defaults to NO FILE LOADED when created)
   * @param parent The (side) panel that holds the button
   */
  private static void initFile(
      JButton button, JTextField text, JPanel parent, GridBagConstraints gbc){

    text.setText("NO FILE LOADED");
    text.setAlignmentX(SwingConstants.CENTER);
    text.setHorizontalAlignment(JTextField.CENTER);

    JScrollPane jsp = new JScrollPane();

    jsp.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
    jsp.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    
    jsp.setViewportView(text);

    gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.NONE;
    parent.add(button, gbc);
    gbc.gridy += 1;
    gbc.fill = GridBagConstraints.HORIZONTAL;
    parent.add(jsp, gbc);
    gbc.gridy += 1;
    
    gbc.fill = GridBagConstraints.NONE;
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

    frame.add( new SensorSuite() );

    frame.pack();
    frame.setVisible(true);
  }

  /**
   * Handles actions when a side-panel button is clicked (file-IO)
   */
  @Override
  public void actionPerformed(ActionEvent e) {

    if ( e.getSource() == clear ) {
      dataBox.clearAllData();
      for (JTextField fn : seedFileNames) {
        fn.setText("NO FILE LOADED");
      }
      for (JTextField fn : respFileNames) {
        fn.setText("NO FILE LOADED");
      }
      return;
    }
    
    for(int i = 0; i < seedLoaders.length; ++i) {
      final int idx = i; // used to play nice with swingworker
      JButton seedButton = seedLoaders[i];
      JButton respButton = respLoaders[i];
      if ( e.getSource() == seedButton ) {
        fc.setCurrentDirectory( new File(seedDirectory) );
        fc.resetChoosableFileFilters();
        fc.setDialogTitle("Load SEED file...");
        int returnVal = fc.showOpenDialog(seedButton);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File file = fc.getSelectedFile();
          seedDirectory = file.getParent();
          
          seedFileNames[i].setText("LOADING: " + file.getName());
          
          // create swingworker to load large files in the background
          SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>(){
            @Override
            public Integer doInBackground() {
              dataBox.setData( idx, file.getAbsolutePath() );
              
              return 0;
            }
            
            public void done() {
              seedFileNames[idx].setText( file.getName() );
              for ( int j = 0; j < tabbedPane.getTabCount(); ++j ) {
                ExperimentPanel ep = 
                    (ExperimentPanel) tabbedPane.getComponentAt(j);
                String[] names = new String[seedFileNames.length];
                for (int k = 0; k < names.length; ++k) {
                  names[k] = seedFileNames[k].getText();
                }
                ep.setDataNames(names);
              }
            }
          };
          // need a new thread so the UI won't lock with big programs
          
          executor.submit(worker);
          return;
        }
      } else if ( e.getSource() == respButton ) {
        // don't need a new thread because resp loading is pretty prompt
        fc.setCurrentDirectory( new File(respDirectory) );
        fc.resetChoosableFileFilters();
        fc.setDialogTitle("Load response file...");
        int returnVal = fc.showOpenDialog(seedButton);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File file = fc.getSelectedFile();
          respDirectory = file.getParent();
          dataBox.setResponse( i, file.getAbsolutePath() );
          respFileNames[i].setText( file.getName() );
        }
        return;
      }
      
    } // end for loop 

    if ( e.getSource() == generate ) {
      this.resetTabPlots();
      return;
    } else if ( e.getSource() == savePDF ) {

      String ext = ".png";
      fc.setCurrentDirectory( new File(saveDirectory) );
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      fc.setDialogTitle("Save plot to PNG...");
      String tStamp = 
          new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format( new Date() );
      fc.setSelectedFile( new File(tStamp+"_ALL.png") );
      int returnVal = fc.showSaveDialog(savePDF);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        saveDirectory = selFile.getParent();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          plotsToPNG(selFile);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
      }
      return;
    }


  }
  
  /**
   * Handles function to create a PNG image with all currently-displayed plots
   * (active experiment and read-in time series data)
   * @param file File (PNG) that image will be saved to
   * @throws IOException
   */
  public void plotsToPNG(File file) throws IOException {
    // using 0s to set image height and width to default values (match window)
    int inHeight = dataBox.getImageHeight();
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

    // for now, it's a png. TODO: change this module write to PDF?

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
