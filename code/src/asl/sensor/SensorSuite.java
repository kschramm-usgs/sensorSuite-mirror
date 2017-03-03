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

import javax.imageio.ImageIO;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.pdmodel.PDPage;
import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.apache.pdfbox.pdmodel.graphics.image.LosslessFactory;
import org.apache.pdfbox.pdmodel.graphics.image.PDImageXObject;

// TODO: move file operations to the datapanel object?

/**
 * Main window of the sensor test program and the program's launcher
 * Mainly used for handling the input (inputpanel) and output (experimentpanel)
 * GUI frames and making sure they fit togther and cooperate
 * @author akearns
 *
 */
public class SensorSuite extends JPanel 
                         implements ActionListener, ChangeListener {

  /**
   * 
   */
  private static final long serialVersionUID = 2866426897343097822L;
 

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
   * Create space with the given dimensions as a buffer between plots
   * when writing the plots to an image file
   * @param width The width of the spacer
   * @param height The height of the spacer
   * @return A blank BufferedImage that can be concatenated with other plots
   */
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

  private JFileChooser fc; // loads in files based on parameter
  
  private InputPanel inputPlots;


  private JTabbedPane tabbedPane; // holds set of experiment panels

  private JButton generate, savePDF; // run all calculations

  // used to store current directory locations
  private String saveDirectory = System.getProperty("user.home");


  /**
   * Creates the main window of the program when called
   * (Three main panels: the top panel for displaying the results
   * of sensor tests; the lower panel for displaying plots of raw data from
   * miniSEED files; the side panel for most file-IO operations
   */
  public SensorSuite() {

    super();
    
    // set up experiment panes in a tabbed pane
    tabbedPane = new JTabbedPane();
    
    for ( ExperimentEnum exp : ExperimentEnum.values() ) {
      JPanel tab = ExperimentPanelFactory.createPanel(exp);
      tabbedPane.addTab( exp.getName(), tab);
    }
    
    Dimension d = tabbedPane.getPreferredSize();
    d.setSize( d.getWidth() * 1.5, d.getHeight() );
    tabbedPane.setPreferredSize(d);
    
    tabbedPane.addChangeListener(this);
    
    inputPlots = new InputPanel();
    inputPlots.addChangeListener(this);
    
    // experiments on left, input on the right; split to allow resizing
    JSplitPane mainSplit = 
        new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, true);
        // boolean allows redrawing immediately on resize
    mainSplit.setLeftComponent(tabbedPane);
    mainSplit.setRightComponent(inputPlots);
    // set the left-pane to resize more when window is horizontally stretched
    mainSplit.setResizeWeight(.5);
    mainSplit.setOneTouchExpandable(true);
    
    // we want to make sure the split pane fills the window
    this.setLayout( new GridBagLayout() );
    GridBagConstraints c = new GridBagConstraints();
    c.gridx = 0; c.gridy = 0;
    c.weightx = 1; c.weighty = 1;
    c.gridwidth = 2;
    c.fill = GridBagConstraints.BOTH;
    
    this.add(mainSplit, c);
    
    c.fill = GridBagConstraints.BOTH;
    c.gridx = 0; c.gridy = 1;
    c.weightx = 1.0; c.weighty = 0.0;
    c.gridwidth = 1;
    
    // now add the buttons
    savePDF = new JButton("Save input and output plots (PDF)");
    savePDF.setEnabled(true); // TODO: change this back?
    savePDF.addActionListener(this);
    this.add(savePDF, c);
    c.gridx += 1;

    generate = new JButton("Generate test result");
    generate.setEnabled(false);
    generate.addActionListener(this);
    d = generate.getPreferredSize();
    d.setSize( d.getWidth(), d.getHeight() * 2 );
    generate.setPreferredSize(d);
    this.add(generate, c);
    
    fc = new JFileChooser();
    
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
    ExperimentEnum em = ep.expType;
    inputPlots.showDataNeeded(em);

  }

  /**
   * Handles actions when a side-panel button is clicked; loading SEED and RESP
   * files and writing the combined set of plots to a file. Because SEED files
   * can be large and take a long time to read in, that operation runs on a
   * separate thread.
   */
  @Override
  public void actionPerformed(ActionEvent e) {


    if ( e.getSource() == generate ) {
      SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {
        
        @Override
        public Integer doInBackground() {
          try{
            resetTabPlots();
          } catch (Exception e) {
            e.printStackTrace();
          }
          return 0;
        }
        
      };
      worker.execute();
      return;
    } else if ( e.getSource() == savePDF ) {

      String ext = ".pdf";
      fc.setCurrentDirectory( new File(saveDirectory) );
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PDF file (.pdf)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      fc.setDialogTitle("Save plot to PNG...");
      String tStamp = 
          new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format( new Date() );
      fc.setSelectedFile( new File(tStamp+"_ALL"+ext) );
      int returnVal = fc.showSaveDialog(savePDF);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        saveDirectory = selFile.getParent();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          plotsToPDF(selFile);
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

    // just write the bufferedimage to file
    ImageIO.write( getCompiledImage(), "png", file );

  }
  
  /**
   * Produces a buffered image of all active charts
   * @return BufferedImage that can be written to file
   */
  public BufferedImage getCompiledImage() {
    // TODO: split this section into its own method (redundant with PNG save)
    int inHeight = inputPlots.getImageHeight() * 2;
    int width = 1280;
    BufferedImage inPlot = inputPlots.getAsImage(width, inHeight);
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();

    BufferedImage outPlot = ep.getAsImage(width, 960);

    // int width = Math.max( inPlot.getWidth(), outPlot.getWidth() );
    // 5px tall buffer used to separate result plot from inputs
    // BufferedImage space = getSpace(width, 0);
    int height = inPlot.getHeight() + outPlot.getHeight();

    // System.out.println(space.getHeight());

    BufferedImage returnedImage = new BufferedImage(width, height, 
        BufferedImage.TYPE_INT_ARGB);

    Graphics2D combined = returnedImage.createGraphics();
    combined.drawImage(outPlot, null, 0, 0);
    combined.drawImage( inPlot, null, 0, 
        outPlot.getHeight() );
    combined.dispose();
    
    return returnedImage;
  }
  
  /**
   * Output the currently displayed plots as a PDF file.
   * @param file Filename to write to
   * @throws IOException If the file cannot be written
   */
  public void plotsToPDF(File file) throws IOException {
    
    BufferedImage toFile = getCompiledImage();
    
    // START OF UNIQUE CODE FOR PDF CREATION HERE
    PDDocument pdf = new PDDocument();
    PDRectangle rec = 
        new PDRectangle( (float) toFile.getWidth(), 
                         (float) toFile.getHeight() );
    PDPage page = new PDPage(rec);
    pdf.addPage(page);
    PDImageXObject  pdImageXObject = 
        LosslessFactory.createFromImage(pdf, toFile);
    PDPageContentStream contentStream = 
        new PDPageContentStream(pdf, page, 
                                PDPageContentStream.AppendMode.OVERWRITE, 
                                true, false);
    contentStream.drawImage( pdImageXObject, 0, 0, 
        toFile.getWidth(), toFile.getHeight() );
    contentStream.close();
    pdf.save( file );
    pdf.close();
  }

  /**
   * Resets plot data and gets the inputted data to send to experiments
   */
  private void resetTabPlots() {
    
    // pass the inputted data to the panels that handle them
    DataStore ds = inputPlots.getData();
    // update the input plots to show the active region being calculated
    inputPlots.showRegionForGeneration();
    
    // now, update the data
    ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
    ep.updateData(ds);
    
    savePDF.setEnabled(true);
  }

  /**
   * Checks when input panel gets new data or active experiment changes
   * to determine whether or not the experiment can be run yet
   */
  @Override
  public void stateChanged(ChangeEvent e) {
    if ( e.getSource() == inputPlots || e.getSource() == tabbedPane ) {
      ExperimentPanel ep = (ExperimentPanel) tabbedPane.getSelectedComponent();
      DataStore ds = inputPlots.getData();
      // TODO: refactor this based around first X inputs
      boolean canGenerate = ep.haveEnoughData(ds);
      ExperimentEnum em = ep.expType;
      inputPlots.showDataNeeded(em);
      generate.setEnabled(canGenerate);
    }
  }


}
