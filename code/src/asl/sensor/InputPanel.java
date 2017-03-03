package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import javax.imageio.ImageIO;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.SwingWorker;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.JTextComponent;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;


/**
 * Panel used to hold the plots for the files taken in as input
 * Handles the UI for loading SEED and RESP files, plotting their data
 * and selecting regions of that data, as well as outputting the plots to PNG
 * @author akearns
 *
 */
public class InputPanel 
extends JPanel 
implements ActionListener, ChangeListener {

  /**
   * auto-generated serialization UID
   */
  private static final long serialVersionUID = -7302813951637543526L;
  
  /**
   * Default height of image produced by the save-as-image function
   * (each chart is 240 pixels tall)
   */
  public static final int IMAGE_HEIGHT = 240;
  public static final int IMAGE_WIDTH = 480;
  
  public static final int FILE_COUNT = DataStore.FILE_COUNT;
  
  private static final int MARGIN = 10; // min space of the two sliders
  private static final int SLIDER_MAX = 10000;
  
  /**
   * Gets the value of start or end time from slider value and DataBlock
   * @param db DataBlock corresponding to one of the plots
   * @param sliderValue Value of starting or ending time slider [0-SLIDER_MAX]
   * @return Long that represents start or end time matching slider's value
   */
  public static long getMarkerLocation(DataBlock db, int sliderValue) {
    long start = db.getStartTime();
    long len = (db.getInterval()) * db.size();
    return start + (sliderValue * len) / SLIDER_MAX;
  }
  private DataStore ds;
  private DataStore zooms;
  private ChartPanel[] chartPanels = new ChartPanel[FILE_COUNT];
  private Color[] defaultColor = {
          ChartColor.LIGHT_RED, 
          ChartColor.LIGHT_BLUE, 
          ChartColor.LIGHT_GREEN };
  private JButton save;
  private JButton zoomIn;
  private JButton zoomOut;
  private JButton clearAll;
  private JFileChooser fc;
  private JPanel allCharts; // parent of the chartpanels, used for image saving
  private JSlider leftSlider;
  private JSlider rightSlider;
      
  private JButton[] seedLoaders  = new JButton[FILE_COUNT];
  private JTextComponent[] seedFileNames = 
      new JTextComponent[FILE_COUNT];
  private JButton[] respLoaders  = new JButton[FILE_COUNT];
  private JTextComponent[] respFileNames = 
      new JTextComponent[FILE_COUNT];
  
  private JButton[] clearButton = new JButton[FILE_COUNT];
  // used to store current directory locations
  private String seedDirectory = "data";
  private String respDirectory = "responses";

  
  private String saveDirectory = System.getProperty("user.home");
  
  private JPanel makeChartSubpanel(int i) {
    
    JPanel chartSubpanel = new JPanel();
    chartSubpanel.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    
    gbc.weightx = 1.0;
    gbc.weighty = 1.0;
    gbc.gridy = 0;
    gbc.anchor = GridBagConstraints.CENTER;
    
    instantiateChart(i);
    
    /*
    chartPanels[i].setMaximumSize(
        new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE) );
    */
    
    chartPanels[i].setMouseZoomable(true);
    
    seedLoaders[i] = new JButton( "Load SEED file " + (i+1) );
    seedLoaders[i].addActionListener(this);
    
    JTextField text = new JTextField( "NO FILE LOADED" );
    text.setHorizontalAlignment(SwingConstants.CENTER);
    seedFileNames[i] = text;
    seedFileNames[i].setEditable(false);
   
    respLoaders[i] = new JButton( "Load RESP file " + (i+1) );
    respLoaders[i].addActionListener(this);
    
    text = new JTextField( "NO FILE LOADED" );
    text.setHorizontalAlignment(SwingConstants.CENTER);
    respFileNames[i] = text;
    respFileNames[i].setEditable(false);
    
    clearButton[i] = new JButton( "Clear data " + (i+1) );
    
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.gridwidth = 1; gbc.gridheight = 5;
    gbc.weightx = 1; gbc.weighty = 1;
    gbc.fill = GridBagConstraints.BOTH;
    chartSubpanel.add(chartPanels[i], gbc);
    
    Dimension d = chartPanels[i].getPreferredSize();
    d.setSize( d.getWidth() / 1.5 , d.getHeight()/1.5 );
    chartPanels[i].setPreferredSize(d);
    
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridx = 1;
    gbc.gridwidth = 1; gbc.gridheight = 1;
    gbc.weightx = 0; gbc.weighty = 0.25;
    chartSubpanel.add(seedLoaders[i], gbc);
    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.weighty += 1;
    gbc.gridy += 1;
    JScrollPane jsp = new JScrollPane(seedFileNames[i]);
    jsp.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
    jsp.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    chartSubpanel.add(jsp, gbc);
    
    gbc.fill = GridBagConstraints.NONE;
    gbc.weighty = 0.25;
    gbc.gridy += 1;
    chartSubpanel.add(respLoaders[i], gbc);

    
    gbc.fill = GridBagConstraints.BOTH;
    gbc.weighty = 1;
    gbc.gridy += 1;
    jsp = new JScrollPane(respFileNames[i]);
    jsp.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
    jsp.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    chartSubpanel.add(jsp, gbc);
    
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.weighty = 0;
    gbc.gridy += 1;
    clearButton[i] = new JButton( "Clear data " + (i+1) );
    clearButton[i].setOpaque(true);
    clearButton[i].setBackground( Color.RED.darker() );
    clearButton[i].addActionListener(this);
    clearButton[i].setEnabled(false);
    chartSubpanel.add(clearButton[i], gbc);
    
    return chartSubpanel;
  }
  
  
  /**
   * Creates a new data panel -- instantiates each chart, to be populated with
   * data when a file is loaded in. Also creates a save button for writing all
   * the inputted data plots into a single PNG file.
   */
  public InputPanel() {
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
   
    ds = new DataStore();
    zooms = new DataStore();
    
    gbc.weightx = 1.0;
    gbc.weighty = 1.0;
    gbc.gridy = 0;
    gbc.gridwidth = 8;
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.fill = GridBagConstraints.BOTH;
    
    JScrollPane jsp = new JScrollPane();
    jsp.setVerticalScrollBarPolicy(
        ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED);
    jsp.setHorizontalScrollBarPolicy(
        ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
    
    Container cont = new Container();
    cont.setLayout( new BoxLayout(cont, BoxLayout.Y_AXIS) );
    Dimension minDim = new Dimension(0, 0);
    
    for (int i = 0; i < FILE_COUNT; ++i) {
      Dimension d = cont.getPreferredSize();
      JPanel chartSubpanel = makeChartSubpanel(i);
      Dimension d2 = chartSubpanel.getPreferredSize();
      minDim = chartSubpanel.getMinimumSize();
      
      d2.setSize( d.getWidth(), d2.getHeight() );
      
      chartSubpanel.setSize( d2 );
      cont.add( makeChartSubpanel(i) );
      // gbc.gridy += 1;
      
    }
    
    jsp.getViewport().setView(cont);
    jsp.setVisible(true);
    jsp.setMinimumSize(minDim);
    
    this.add(jsp, gbc);
    gbc.gridy += 1;
    
    // set size so that the result pane isn't distorted on window launch
    Dimension d = jsp.getPreferredSize();
    d.setSize( d.getWidth() + 5, d.getHeight()* 2 );
    this.setPreferredSize(d);
    
    leftSlider = new JSlider(0, SLIDER_MAX, 0);
    leftSlider.setEnabled(false);
    leftSlider.addChangeListener(this);
    
    rightSlider = new JSlider(0, SLIDER_MAX, SLIDER_MAX);
    rightSlider.setEnabled(false);
    rightSlider.addChangeListener(this);
    
    zoomIn = new JButton("Zoom in (on selection)");
    zoomIn.addActionListener(this);
    zoomIn.setEnabled(false);
    
    zoomOut = new JButton("Zoom out (show all)");
    zoomOut.addActionListener(this);
    zoomOut.setEnabled(false);
    
    save = new JButton("Save input (PNG)");
    save.addActionListener(this);
    save.setEnabled(false);
    
    clearAll = new JButton("Clear ALL data");
    clearAll.setOpaque(true);
    clearAll.setBackground( Color.RED.darker() );
    clearAll.addActionListener(this);
    clearAll.setEnabled(false);
    
    int yOfClear = gbc.gridy;
    
    //this.add(allCharts);
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridy += 2;
    gbc.weighty = 0;
    gbc.weightx = 1;
    gbc.gridwidth = 1;
    gbc.gridx = 0;

    this.add(leftSlider, gbc);
    
    gbc.gridx += 1;
    this.add(rightSlider, gbc);
    
    // now we can add the space between the last plot and the save button
    //this.add( Box.createVerticalStrut(5) );
    
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1;
    gbc.weighty = 0; gbc.weightx = 0;
    
    gbc.gridx = 0;
    gbc.anchor = GridBagConstraints.EAST;
    this.add(zoomIn, gbc);
    
    gbc.gridx += 1;
    gbc.anchor = GridBagConstraints.WEST;
    this.add(zoomOut, gbc);

    
    gbc.gridwidth = 1;
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.gridx = 7;
    gbc.gridy = yOfClear+1;
    gbc.gridheight = 2;
    gbc.fill = GridBagConstraints.BOTH;

    this.add(save, gbc);

    
    gbc.gridy += 2;
    gbc.gridheight = GridBagConstraints.REMAINDER;

    this.add(clearAll, gbc);

    
    //this.add(buttons);
    
    fc = new JFileChooser();
    
  }
  
  /**
   * Dispatches commands when interface buttons are clicked.
   * When the save button is clicked, dispatches the command to save plots as
   * an image. When the zoom buttons are clicked, scales the plot to only
   * contain the data within a specific range. Resets plots and removes 
   * underlying data when the clear buttons are clicked. Prompts user to
   * load in a file for miniseed and resp data when the corresponding
   * buttons are clicked; because seed files can be large and contian a lot
   * of data to plot, runs seed-loading code backend in a separate thread.
   * 
   * When new data is loaded in, this also fires a change event to any listeners
   */
  @Override
  public void actionPerformed(ActionEvent e) {

    for (int i = 0; i < FILE_COUNT; ++i) {
      JButton clear = clearButton[i];
      JButton seed = seedLoaders[i];
      JButton resp = respLoaders[i];
      
      if( e.getSource() == clear ) {
        instantiateChart(i);
        ds.removeData(i);
        clear.setEnabled(false);
        seedFileNames[i].setText("NO FILE LOADED");
        respFileNames[i].setText("NO FILE LOADED");
        
        // plot all valid range of loaded-in data or else
        // disable clearAll button if there's no other data loaded in
        clearAll.setEnabled( ds.isAnythingSet() );
        zooms = new DataStore(ds);
        zooms.trimToCommonTime();
        
        showRegionForGeneration();
        
        fireStateChanged();
      }
      
      if ( e.getSource() == seed ) {
        loadData(i, seed);
      }
      
      if ( e.getSource() == resp ) {
        // don't need a new thread because resp loading is pretty prompt
        Set<String> respFilenames = new HashSet<String>();
        ClassLoader cl = InputPanel.class.getClassLoader();
        
        // there's no elegant way to extract responses other than to
        // load in their names from a list and then grab them as available
        // correspondingly, this means adding response files to this program
        // requires us to add their names to this file
        // There may be other possibilities but they are more complex and
        // tend not to work the same way between IDE and launching a jar
        
        InputStream respRead = cl.getResourceAsStream("responses.txt");
        BufferedReader respBuff = 
            new BufferedReader( new InputStreamReader(respRead) );

        try {
          String name;
          name = respBuff.readLine();
          while (name != null) {
            respFilenames.add(name);
            name = respBuff.readLine();
          }
          respBuff.close();
        } catch (IOException e2) {
          // TODO Auto-generated catch block
          e2.printStackTrace();
        }
        
        List<String> names = 
            new ArrayList<String>(respFilenames);
        
        String custom = "Load custom response...";
        
        names.add(custom);
        
        JDialog dialog = new JDialog();
        Object result = JOptionPane.showInputDialog(
            dialog,
            "Select a response to load:",
            "RESP File Selection",
            JOptionPane.PLAIN_MESSAGE,
            null, names.toArray(),
            names.get( names.size() - 1 ) );
        
        String resultStr = (String) result;
        
        if (resultStr == null) {
          return;
        }
        
        if ( respFilenames.contains(resultStr) ) {
          // TODO: load response mappings (need responses first)
          final String fname = resultStr;
          InputStream is = cl.getResourceAsStream(fname);
          BufferedReader fr = new BufferedReader( new InputStreamReader(is) );
          try {
            InstrumentResponse ir = new InstrumentResponse(fr);
            ds.setResponse(i, ir);
            zooms.setResponse(i, ir);
            respFileNames[i].setText( resultStr );
            clear.setEnabled(true);
            clearAll.setEnabled(true);
            
            fireStateChanged();
          } catch (IOException e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
          }
        } else {
          fc.setCurrentDirectory( new File(respDirectory) );
          fc.resetChoosableFileFilters();
          fc.setDialogTitle("Load response file...");
          int returnVal = fc.showOpenDialog(resp);
          if (returnVal == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            respDirectory = file.getParent();
            ds.setResponse(i, file.getAbsolutePath() );
            zooms.setResponse(i, file.getAbsolutePath() );
            respFileNames[i].setText( file.getName() );
            clear.setEnabled(true);
            clearAll.setEnabled(true);
            
            fireStateChanged();
          }
        }
        

        return;
      }
      
    }
    
    if ( e.getSource() == clearAll) {
      clearAllData();
      
      fireStateChanged();
      
      return;
    }
    
    
    if ( e.getSource() == save ) {
      String ext = ".png";
      fc = new JFileChooser();
      fc.setCurrentDirectory( new File(saveDirectory) );
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      int returnVal = fc.showSaveDialog(save);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        saveDirectory = selFile.getParent();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          BufferedImage bi = getAsImage( 640, getImageHeight() );
          
          ImageIO.write(bi,"png",selFile);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
        
      }
      return;
    }
    
    if ( e.getSource() == zoomIn ) {
      
      showRegionForGeneration();
      return;
    }
    
    if ( e.getSource() == zoomOut ) {
      
      // restore original loaded datastore
      zooms = new DataStore(ds);
      zooms.trimToCommonTime();
      for (int i = 0; i < FILE_COUNT; ++i) {
        if ( !zooms.blockIsSet(i) ) {
          continue;
        }
        resetPlotZoom(i);
      }
      
      leftSlider.setValue(0); rightSlider.setValue(SLIDER_MAX);
      setVerticalBars();
      zoomOut.setEnabled(false);
      return;
    }
  }
  
  /**
   * Used to add objects to the list that will be informed when data is loaded
   * @param listener
   */
  public void addChangeListener(ChangeListener listener) {
    listenerList.add(ChangeListener.class, listener);
  }
  
  /**
   * Resets the data and blanks out all charts
   */
  public void clearAllData() {
    ds = new DataStore();
    zooms = new DataStore();
    
    zoomIn.setEnabled(false);
    zoomOut.setEnabled(false);
    
    leftSlider.setEnabled(false);
    rightSlider.setEnabled(false);
    
    clearAll.setEnabled(false);
    save.setEnabled(false);
    
    for (JTextComponent fn : seedFileNames) {
      fn.setText("NO FILE LOADED");
    }
    for (JTextComponent fn : respFileNames) {
      fn.setText("NO FILE LOADED");
    }
    
    for (int i = 0; i < chartPanels.length; ++i) {
      clearButton[i].setEnabled(false);
      instantiateChart(i);
    }
    
  }
  
  /**
   * Informs listening objects that the state of the inputs has changed
   * This is done when new seed or resp data has been loaded in, mainly
   * to tell whether enough data exists to run one of the experiments
   */
  protected void fireStateChanged() {
    ChangeListener[] lsners = listenerList.getListeners(ChangeListener.class);
    if (lsners != null && lsners.length > 0) {
      ChangeEvent evt = new ChangeEvent(this);
      for (ChangeListener lsnr : lsners) {
        lsnr.stateChanged(evt);
      }
    }
  }
  
  /**
   * Return this panel's charts as a single buffered image
   * @return Buffered image of the plots, writeable to file
   */
  public BufferedImage getAsImage() {
    // if all 3 plots are set, height of panel is height of image
    int height = getImageHeight();
    // otherwise, we only use the height of images actually set
    return getAsImage( IMAGE_WIDTH, height );
  }
  
  /**
   * Return this panel's charts as a single buffered image 
   * with specified dimensions
   * @param width Width of returned image
   * @param height Height of returned image
   * @return Buffered image of the plots, writeable to file
   */
  public BufferedImage getAsImage(int width, int height) {
    
    // int shownHeight = allCharts.getHeight();
    
    // width = Math.min( width, chartPanels[0].getWidth() );
    // height = Math.max( height, shownHeight );
    
    int loaded = ds.numberOfBlocksSet();
    // cheap way to make sure height is a multiple of the chart count
    height = (height*loaded)/loaded;
    
    int chartHeight = height/loaded;
    
    Dimension outSize = new Dimension(width, height);

    JPanel toDraw = new JPanel(); // what we're going to draw to image
    toDraw.setSize(outSize);
    toDraw.setPreferredSize(outSize);
    toDraw.setMinimumSize(outSize);
    toDraw.setMaximumSize(outSize);
    
    toDraw.setLayout( new BoxLayout(toDraw, BoxLayout.Y_AXIS) );

    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( !ds.blockIsSet(i) ) {
        continue;
      }
      ChartPanel cp = chartPanels[i];
      // toDraw.add( Box.createVerticalStrut(5) );
      Dimension chartSize = new Dimension(width, chartHeight);
      ChartPanel outPanel = new ChartPanel( cp.getChart() );
      outPanel.setSize(chartSize);
      outPanel.setPreferredSize(chartSize);
      outPanel.setMinimumSize(chartSize);
      outPanel.setMaximumSize(chartSize);
      outPanel.setMinimumDrawHeight(chartHeight);
      toDraw.add(outPanel);
    }
    
    // used to make sure that everything is laid out correctly when we save
    // (forces the Java windowing tools to respect the specified layout above)
    // we do this since the panel is instantiated here, not displayed, and
    // is built from multiple subcomponents, unlike experimentpanel
    // Before the frame exists, Java tends to ignore any layout instructions
    JFrame jw = new JFrame();
    jw.add(toDraw);
    jw.pack();
    
    BufferedImage bi = new BufferedImage(
        width, 
        height, 
        BufferedImage.TYPE_INT_ARGB);

    Graphics2D g = bi.createGraphics();
    toDraw.printAll(g);
    g.dispose();
    
    return bi;
  }
  
  
  /**
   * Returns the selected region of underlying DataStore, to be fed 
   * into experiments for processing (the results of which will be plotted)
   * When this function is called, the graphs zoom to the currently active
   * range to display selected by the sliders, which is also the range
   * passed into an experiment
   * @return A DataStore object (contains arrays of DataBlocks & Responses)
   */
  public DataStore getData() {
    
    if ( ds.numberOfBlocksSet() < 1 ) {
      return new DataStore();
    }
    

    showRegionForGeneration();
    zooms.matchIntervals();
    zooms.trimToCommonTime();
    return zooms;
  }

  /**
   * Gets the height of resulting image of plots given default parameters,
   * so that it only needs to fit the plots that have data in them 
   * @return height of image to output
   */
  public int getImageHeight() {
    return IMAGE_HEIGHT * ds.numberOfBlocksSet();
  }
  
  /**
   * Returns a default image width for writing plots to file
   * @return width of image to output
   */
  public int getImageWidth()  {
    return allCharts.getWidth();
  }
  
  /**
   * Instantiates the underlying chart of a chartpanel with default data
   * @param idx Index of the chartpanel to instantiate
   */
  private void instantiateChart(int idx) {
    JFreeChart chart = ChartFactory.createXYLineChart(
        "SEED input " + (idx + 1),
        "Time",
        "Counts",
        new XYSeriesCollection(),
        PlotOrientation.VERTICAL,
        false, false, false);
    
    if (chartPanels[idx] == null) {
      chartPanels[idx] = new ChartPanel(chart);
    } else {
      chartPanels[idx].setChart(chart);

    }
    chartPanels[idx].setMouseZoomable(true);
  }
  
  /**
   * Used to remove an object from the list of those informed when
   * data is loaded in or cleared out
   * @param listener
   */
  public void removeChangeListener(ChangeListener listener) {
    listenerList.remove(ChangeListener.class, listener);
  }
  
  /**
   * Does the work to reset the zoom of a chart when the zoom button is hit
   * @param idx Index of appropriate chart/panel
   */
  private void resetPlotZoom(int idx) {
    XYPlot xyp = chartPanels[idx].getChart().getXYPlot();
    XYSeriesCollection xys = new XYSeriesCollection();
    xys.addSeries( zooms.getBlock(idx).toXYSeries() );
    xyp.setDataset( xys );
    xyp.getRenderer().setSeriesPaint(0,
        defaultColor[idx % defaultColor.length]);
    if ( xyp.getSeriesCount() > 1 ) {
      throw new RuntimeException("TOO MUCH DATA");
    }
    chartPanels[idx].repaint();
  }
  
  /**
   * Parent function to load in a response file to this panel's
   * DataStore object
   * @param idx Index of the file/chart this response corresponds to
   * @param filepath Full address of the file to be loaded in
   */
  public void setResponse(int idx, String filepath) {
    
    ds.setResponse(idx, filepath);
    zooms.setResponse(idx, filepath);
  }
  
  /**
   * Displays the range set by the sliders using
   * vertical bars at the min and max values
   */
  public void setVerticalBars() {
    
    // zooms.trimToCommonTime();
    
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( !ds.blockIsSet(i) ) {
        continue;
      }
      
      int leftValue = leftSlider.getValue();
      int rightValue = rightSlider.getValue();
      
      XYPlot xyp = (XYPlot) chartPanels[i].getChart().getPlot();
      xyp.clearDomainMarkers();
      
      DataBlock db = zooms.getBlock(i);
      
      long startMarkerLocation = getMarkerLocation(db, leftValue);
      long endMarkerLocation = getMarkerLocation(db, rightValue);
      
      // divide by 1000 here to get time value in ms
      Marker startMarker = new ValueMarker(startMarkerLocation/1000);
      startMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      Marker endMarker = new ValueMarker(endMarkerLocation/1000);
      endMarker.setStroke( new BasicStroke( (float) 1.5 ) );
      
      xyp.addDomainMarker(startMarker);
      xyp.addDomainMarker(endMarker);
      
      chartPanels[i].repaint();
    }
    
    
  }

  /**
   * Zooms in on the current range of data, which will be passed into
   * backend functions for experiment calculations
   */
  public void showRegionForGeneration() {
    
    if ( ds.numberOfBlocksSet() < 1 ) {
      return;
    }
    
    // get (any) loaded data block to map slider to domain boundary
    // all data should have the same range
    DataBlock db = zooms.getXthLoadedBlock(1);

    if ( leftSlider.getValue() != 0 || rightSlider.getValue() != SLIDER_MAX ) {
      long start = getMarkerLocation(db, leftSlider.getValue() );
      long end = getMarkerLocation(db, rightSlider.getValue() );
      zooms = new DataStore(ds, start, end);
      leftSlider.setValue(0); rightSlider.setValue(SLIDER_MAX);
      zoomOut.setEnabled(true);
    }
    
    for (int i = 0; i < FILE_COUNT; ++i) {
      if ( !zooms.blockIsSet(i) ) {
        continue;
      }
      resetPlotZoom(i);
    }

    setVerticalBars();
    
  }

  @Override
  /**
   * Handles changes in value by the sliders below the charts
   */
  public void stateChanged(ChangeEvent e) {
    
    int leftSliderValue = leftSlider.getValue();
    int rightSliderValue = rightSlider.getValue();
    
    // probably can refactor this
    // the conditionals are effectively the same
    
    if ( e.getSource() == leftSlider ) {
      if (leftSliderValue > rightSliderValue || 
          leftSliderValue + MARGIN > rightSliderValue) {
        leftSliderValue = rightSliderValue - MARGIN;
        if (leftSliderValue < 0) {
          leftSliderValue = 0;
          rightSliderValue = MARGIN;
        }
      }
    } else if ( e.getSource() == rightSlider ) {
      if (rightSliderValue < leftSliderValue ||
          rightSliderValue - MARGIN < leftSliderValue) {
        rightSliderValue = leftSliderValue + MARGIN;
        if (rightSliderValue > SLIDER_MAX) {
          rightSliderValue = SLIDER_MAX;
          leftSliderValue = SLIDER_MAX - MARGIN;
        }
      }
    }
    
    leftSlider.setValue(leftSliderValue);
    rightSlider.setValue(rightSliderValue);
    
    setVerticalBars();
    
  }
  
  private void loadData(final int idx, final JButton seed) {

    fc.setCurrentDirectory( new File(seedDirectory) );
    fc.resetChoosableFileFilters();
    fc.setDialogTitle("Load SEED file...");
    int returnVal = fc.showOpenDialog(seed);
    if (returnVal == JFileChooser.APPROVE_OPTION) {
      final File file = fc.getSelectedFile();
      seedDirectory = file.getParent();
      String oldName = seedFileNames[idx].getText();

      seedFileNames[idx].setText("LOADING: " + file.getName());
      final String filePath = file.getAbsolutePath();
      String filterName = "";
      try {
        Set<String> nameSet = TimeSeriesUtils.getMplexNameSet(filePath);

        if (nameSet.size() > 1) {
          // more than one series in the file? prompt user for it
          String[] names = nameSet.toArray(new String[0]);
          Arrays.sort(names);
          JDialog dialog = new JDialog();
          Object result = JOptionPane.showInputDialog(
              dialog,
              "Select the subseries to load:",
              "Multiplexed File Selection",
              JOptionPane.PLAIN_MESSAGE,
              null, names,
              names[0]);
          if (result instanceof String) {
            filterName = (String) result;
          } else {
            // if the user cnacelled selecting a subseries
            seedFileNames[idx].setText(oldName);
            return;
          }
        } else {
          // just get the first one; it's the only one in the list
          filterName = new ArrayList<String>(nameSet).get(0);
        }

      } catch (FileNotFoundException e1) {
        e1.printStackTrace();
        return;
      }

      final String immutableFilter = filterName;

      // create swingworker to load large files in the background
      SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>(){

        JFreeChart chart;
        boolean caughtException = false;

        @Override
        public Integer doInBackground() {

          try {
            ds.setData(idx, filePath, immutableFilter);
          } catch (RuntimeException e) {
            caughtException = true;
            return 1;
          }

          XYSeries ts = ds.getPlotSeries(idx);
          chart = ChartFactory.createXYLineChart(
              ts.getKey().toString(),
              "Time",
              "Counts",
              new XYSeriesCollection(ts),
              PlotOrientation.VERTICAL,
              false, false, false);

          XYPlot xyp = (XYPlot) chart.getPlot();
          DateAxis da = new DateAxis();
          SimpleDateFormat sdf = new SimpleDateFormat("Y.DDD.HH:mm");
          sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
          da.setLabel("UTC Time (Year.Day.Hour:Minute)");
          Font bold = da.getLabelFont();
          bold = bold.deriveFont(Font.BOLD);
          da.setLabelFont(bold);
          da.setDateFormatOverride(sdf);
          xyp.setDomainAxis(da);
          xyp.getRenderer().setSeriesPaint(0, defaultColor[idx]);

          zooms = new DataStore(ds);
          zooms.matchIntervals();
          zooms.trimToCommonTime();

          return 0;
          // setData(idx, filePath, immutableFilter);
          // return 0;
        }

        public void done() {

          if (caughtException) {
            instantiateChart(idx);
            XYPlot xyp = (XYPlot) chartPanels[idx].getChart().getPlot();
            TextTitle result = new TextTitle();
            String errMsg = "COULD NOT LOAD IN " + file.getName();
            errMsg += "\nTIME RANGE DOES NOT INTERSECT";
            result.setText(errMsg);
            result.setBackgroundPaint(Color.red);
            result.setPaint(Color.white);
            XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
                RectangleAnchor.CENTER);
            xyp.clearAnnotations();
            xyp.addAnnotation(xyt);
            seedFileNames[idx].setText("NO FILE LOADED");
            clearButton[idx].setEnabled(true);
            return;
          }

          seedFileNames[idx].setText("PLOTTING: " + file.getName());

          chartPanels[idx].setChart(chart);
          chartPanels[idx].setMouseZoomable(true);

          clearButton[idx].setEnabled(true);

          for (int i = 0; i < FILE_COUNT; ++i) {
            if ( !zooms.blockIsSet(i) ) {
              continue;
            }
            resetPlotZoom(i);
          }
          
          leftSlider.setValue(0); rightSlider.setValue(SLIDER_MAX);
          setVerticalBars();

          zoomOut.setEnabled(false);
          zoomIn.setEnabled(true);
          leftSlider.setEnabled(true);
          rightSlider.setEnabled(true);
          save.setEnabled(true);
          clearAll.setEnabled(true);

          leftSlider.setValue(0);
          rightSlider.setValue(SLIDER_MAX);
          setVerticalBars();

          seedFileNames[idx].setText( 
              file.getName() + ": " + immutableFilter);

          fireStateChanged();
        }
      };

      worker.execute();
      return;
    }
  }
  
}
