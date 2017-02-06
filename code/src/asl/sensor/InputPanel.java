package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Set;
import java.util.TimeZone;

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
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


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
  
  private DataStore ds;
  private DataStore zooms;
  private ChartPanel[] chartPanels = new ChartPanel[DataStore.FILE_COUNT];
  private boolean[] set = new boolean[DataStore.FILE_COUNT];
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
  private final int MARGIN = 10; // min space of the two sliders
  
  private JButton[] seedLoaders  = new JButton[DataStore.FILE_COUNT];
  private JTextComponent[] seedFileNames = 
      new JTextComponent[DataStore.FILE_COUNT];
  private JButton[] respLoaders  = new JButton[DataStore.FILE_COUNT];
  private JTextComponent[] respFileNames = 
      new JTextComponent[DataStore.FILE_COUNT];
  private JButton[] clearButton = new JButton[DataStore.FILE_COUNT];
  
  // used to store current directory locations
  private String seedDirectory = "data";
  private String respDirectory = "responses";
  private String saveDirectory = System.getProperty("user.home");

  
  
  /**
   * Creates a new data panel -- instantiates each chart, to be populated with
   * data when a file is loaded in. Also creates a save button for writing all
   * the inputted data plots into a single PNG file.
   */
  public InputPanel() {
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
   
    ds = new DataStore();
    zooms = ds;
    
    gbc.weightx = 1.0;
    gbc.weighty = 1.0;
    gbc.gridy = 0;
    gbc.anchor = GridBagConstraints.CENTER;
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      
      gbc.gridy = i * 5;
      
      JFreeChart chart = ChartFactory.createXYLineChart(
          "",
          "",
          "",
          new XYSeriesCollection(),
          PlotOrientation.VERTICAL,
          false, false, false);
      
      chartPanels[i] = new ChartPanel(chart);
      /*
      chartPanels[i].setMaximumSize(
          new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE) );
      */
      
      chartPanels[i].setMouseZoomable(true);
      
      seedLoaders[i] = new JButton( "Load SEED file " + (i+1) );
      
      JTextField text = new JTextField( "NO FILE LOADED" );
      text.setHorizontalAlignment(SwingConstants.CENTER);
      seedFileNames[i] = text;
      seedFileNames[i].setEditable(false);
     
      respLoaders[i] = new JButton( "Load RESP file " + (i+1) );
      
      text = new JTextField( "NO FILE LOADED" );
      text.setHorizontalAlignment(SwingConstants.CENTER);
      respFileNames[i] = text;
      respFileNames[i].setEditable(false);
      
      clearButton[i] = new JButton( "Clear data " + (i+1) );
      
      gbc.gridx = 0;
      gbc.gridwidth = 6; gbc.gridheight = 5;
      gbc.weightx = 1; gbc.weighty = 1;
      gbc.fill = GridBagConstraints.BOTH;
      this.add(chartPanels[i], gbc);
      
      gbc.fill = GridBagConstraints.BOTH;
      gbc.gridx = 7; gbc.gridy = i * 5;
      gbc.gridwidth = 1; gbc.gridheight = 1;
      gbc.weightx = 0; gbc.weighty = 0.25;
      this.add(seedLoaders[i], gbc);
      seedLoaders[i].addActionListener(this);
      
      gbc.fill = GridBagConstraints.BOTH;
      gbc.weighty = 1;
      gbc.gridy += 1;
      JScrollPane jsp = new JScrollPane(seedFileNames[i]);
      jsp.setVerticalScrollBarPolicy(
          ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
      jsp.setHorizontalScrollBarPolicy(
          ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
      this.add(jsp, gbc);
      
      gbc.fill = GridBagConstraints.BOTH;
      gbc.weighty = 0.25;
      gbc.gridy += 1;
      this.add(respLoaders[i], gbc);
      respLoaders[i].addActionListener(this);
      
      gbc.fill = GridBagConstraints.BOTH;
      gbc.weighty = 1;
      gbc.gridy += 1;
      jsp = new JScrollPane(respFileNames[i]);
      jsp.setVerticalScrollBarPolicy(
          ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
      jsp.setHorizontalScrollBarPolicy(
          ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
      this.add(jsp, gbc);
      
      gbc.fill = GridBagConstraints.HORIZONTAL;
      gbc.weighty = 0;
      gbc.gridy += 1;
      clearButton[i] = new JButton("Clear this data");
      clearButton[i].setOpaque(true);
      clearButton[i].setBackground( Color.RED.darker() );
      this.add(clearButton[i], gbc);
      clearButton[i].addActionListener(this);
      clearButton[i].setEnabled(false);
      
    }
    
    // set size so that the result pane isn't distorted on window launch
    Dimension d = chartPanels[0].getPreferredSize();
    int height = (int) d.getHeight();
    int width = (int) d.getWidth();
    
    this.setPreferredSize( new Dimension(width, height*2) );
    
    int yOfClear = gbc.gridy;
    
    //this.add(allCharts);
    gbc.fill = GridBagConstraints.HORIZONTAL;
    gbc.gridy += 2;
    gbc.weighty = 0;
    gbc.weightx = 1;
    gbc.gridwidth = 3;
    gbc.gridx = 0;
    leftSlider = new JSlider(0, 1000, 0);
    leftSlider.setEnabled(false);
    leftSlider.addChangeListener(this);
    this.add(leftSlider, gbc);
    
    gbc.gridwidth = 1;
    gbc.gridx += 1;
    this.add(new JPanel(), gbc);
    
    gbc.gridx += 3;
    rightSlider = new JSlider(0, 1000, 1000);
    rightSlider.setEnabled(false);
    rightSlider.addChangeListener(this);
    this.add(rightSlider, gbc);
    
    // now we can add the space between the last plot and the save button
    //this.add( Box.createVerticalStrut(5) );
    
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1;
    gbc.weighty = 0; gbc.weightx = 0;
    
    gbc.gridwidth = 3;
    
    gbc.gridx = 0;
    gbc.anchor = GridBagConstraints.EAST;
    zoomIn = new JButton("Zoom in (on selection)");
    this.add(zoomIn, gbc);
    zoomIn.addActionListener(this);
    zoomIn.setEnabled(false);

    
    gbc.gridx += 3;
    gbc.anchor = GridBagConstraints.WEST;
    zoomOut = new JButton("Zoom out (show all)");
    this.add(zoomOut, gbc);
    zoomOut.addActionListener(this);
    zoomOut.setEnabled(false);
    
    gbc.gridwidth = 1;
    gbc.anchor = GridBagConstraints.CENTER;
    gbc.gridx = 7;
    gbc.gridy = yOfClear+1;
    gbc.gridheight = 2;
    gbc.fill = GridBagConstraints.BOTH;
    save = new JButton("Save input (PNG)");
    this.add(save, gbc);
    save.addActionListener(this);
    save.setEnabled(false);
    
    gbc.gridy += 2;
    gbc.gridheight = GridBagConstraints.REMAINDER;
    clearAll = new JButton("Clear all data");
    this.add(clearAll, gbc);
    clearAll.setOpaque(true);
    clearAll.setBackground( Color.RED.darker() );
    clearAll.addActionListener(this);
    clearAll.setEnabled(false);
    
    //this.add(buttons);
    
    fc = new JFileChooser();
    
  }
  
  /**
   * Displays the range set by the sliders using
   * vertical bars at the min and max values
   */
  public void setVerticalBars() {
    
    // zooms.trimToCommonTime();
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      if ( null == ds.getBlock(i) ) {
        continue;
      }
      
      int leftValue = leftSlider.getValue();
      int rightValue = rightSlider.getValue();
      
      XYPlot xyp = (XYPlot) chartPanels[i].getChart().getPlot();
      xyp.clearDomainMarkers();
      
      DataBlock db = zooms.getBlock(i);
      
      long startMarkerLocation = getMarkerLocation(db, leftValue);
      long endMarkerLocation = getMarkerLocation(db, rightValue);
      
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
   * Gets the value of start or end time from slider value and DataBlock
   * @param db DataBlock corresponding to one of the plots
   * @param sliderValue Value of starting or ending time slider [0-1000]
   * @return Long that represents start or end time matching slider's value
   */
  public static long getMarkerLocation(DataBlock db, int sliderValue) {
    long start = db.getStartTime();
    long len = (db.getInterval()) * db.size();
    return start + (sliderValue * len) / 1000;
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
   * Returns the selected region of underlying DataStore, to be fed 
   * into experiments for processing (the results of which will be plotted)
   * @return A DataStore object (contains arrays of DataBlocks & Responses)
   */
  public DataStore getData() {
    
    if ( ds.numberOfBlocksLoaded() < 1 ) {
      return new DataStore();
    }
    

    showRegionForGeneration();
    return zooms;
  }

  /**
   * Dispatches commands based on save and zoom buttons clicked
   */
  @Override
  public void actionPerformed(ActionEvent e) {

    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      JButton clear = clearButton[i];
      JButton seed = seedLoaders[i];
      JButton resp = respLoaders[i];
      
      if( e.getSource() == clear ) {
        JFreeChart chart = ChartFactory.createXYLineChart(
            "",
            "",
            "",
            new XYSeriesCollection(),
            PlotOrientation.VERTICAL,
            false, false, false);
        chartPanels[i].setChart(chart);
        ds.removeData(i);
        clear.setEnabled(false);
        seedFileNames[i].setText("NO FILE LOADED");
        respFileNames[i].setText("NO FILE LOADED");
      }
      
      if ( e.getSource() == seed ) {
        final int idx = i;
        fc.setCurrentDirectory( new File(seedDirectory) );
        fc.resetChoosableFileFilters();
        fc.setDialogTitle("Load SEED file...");
        int returnVal = fc.showOpenDialog(seed);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
          File file = fc.getSelectedFile();
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

          final File immutableFile = file;
          final String immutableFilter = filterName;
          
          // create swingworker to load large files in the background
          SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>(){
            
            JFreeChart chart;
            
            @Override
            public Integer doInBackground() {
              // generate.setEnabled(false);
              ds.setData(idx, filePath, immutableFilter);
              XYSeries ts = ds.getBlock(idx).toXYSeries();
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

              return 0;
              // setData(idx, filePath, immutableFilter);
              // return 0;
            }
            
            public void done() {
              
              zooms = ds;
              for (int j = 0; j < DataStore.FILE_COUNT; ++j) {
                if (zooms.getBlock(j) == null) {
                  continue;
                }
                resetPlotZoom(j);
              }
              
              if (chart == null) {
                set[idx] = false;
                return;
              }
              
              chartPanels[idx].setChart(chart);
              chartPanels[idx].setMouseZoomable(true);

              clearButton[idx].setEnabled(true);
              
              showRegionForGeneration();
              
              set[idx] = true;

              zoomIn.setEnabled(true);
              leftSlider.setEnabled(true);
              rightSlider.setEnabled(true);
              save.setEnabled(true);
              clearAll.setEnabled(true);
              
              leftSlider.setValue(0);
              rightSlider.setValue(1000);
              setVerticalBars();
              
              seedFileNames[idx].setText( 
                  immutableFile.getName() + ": " + immutableFilter);
              // generate.setEnabled(true); TODO callback somehow
            }
          };
          // need a new thread so the UI won't lock with big programs
          
          new Thread(worker).run();
          return;
        }
      }
      
      if ( e.getSource() == resp ) {
        // don't need a new thread because resp loading is pretty prompt
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
          clearAll.setEnabled(true);
        }
        return;
      }
      
    }
    
    if ( e.getSource() == clearAll) {
      clearAllData();
      
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
      zooms = ds;
      for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
        if (zooms.getBlock(i) == null) {
          continue;
        }
        this.resetPlotZoom(i);
      }
      
      leftSlider.setValue(0); rightSlider.setValue(1000);
      setVerticalBars();
      zoomOut.setEnabled(false);
      return;
    }
  }
  
  public void showRegionForGeneration() {
    
    if ( ds.numberOfBlocksLoaded() < 1 ) {
      return;
    }
    
    // get (any) loaded data block to map slider to domain boundary
    // all data should have the same range
    DataBlock db = zooms.getXthLoadedBlock(1);

    if ( leftSlider.getValue() != 0 || rightSlider.getValue() != 1000 ) {
      long start = getMarkerLocation(db, leftSlider.getValue() );
      long end = getMarkerLocation(db, rightSlider.getValue() );
      zooms = new DataStore(ds, start, end);
      for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
        if (zooms.getBlock(i) == null) {
          continue;
        }
        resetPlotZoom(i);
      }
      leftSlider.setValue(0); rightSlider.setValue(1000);
    }
    

    setVerticalBars();
    zoomOut.setEnabled(true);
    
  }
  
  /**
   * Does the work to reset the zoom of each chart
   * @param idx Index of appropriate chart/panel
   */
  private void resetPlotZoom(int idx) {
    XYPlot xyp = chartPanels[idx].getChart().getXYPlot();
    XYSeriesCollection xys = new XYSeriesCollection();
    xys.addSeries( zooms.getBlock(idx).toXYSeries() );
    xyp.setDataset( xys );
    xyp.getRenderer().setSeriesPaint(0, defaultColor[idx]);
    if ( xyp.getSeriesCount() > 1 ) {
      throw new RuntimeException("TOO MUCH DATA");
    }
    chartPanels[idx].repaint();
  }
  
  public int getImageHeight() {
    return IMAGE_HEIGHT * ds.numberOfBlocksLoaded();
  }
  
  public int getImageWidth()  {
    return allCharts.getWidth();
  }
  
  /**
   * Return this panel's charts as a single buffered image
   * @return Buffered image of the plots, writeable to file
   */
  public BufferedImage getAsImage() {
    // if all 3 plots are set, height of panel is height of image
    int height = getImageHeight();
    // otherwise, we only use the height of images actually set
    return getAsImage( allCharts.getWidth(), height );
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
    
    width = Math.min( width, chartPanels[0].getWidth() );
    // height = Math.max( height, shownHeight );
    
    int loaded = ds.numberOfBlocksLoaded();
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

    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
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

  @Override
  /**
   * Handles changes in value by the sliders below the charts
   */
  public void stateChanged(ChangeEvent e) {
    // TODO Auto-generated method stub
    
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
        if (rightSliderValue > 1000) {
          rightSliderValue = 1000;
          leftSliderValue = 1000-MARGIN;
        }
      }
    }
    
    leftSlider.setValue(leftSliderValue);
    rightSlider.setValue(rightSliderValue);
    
    setVerticalBars();
    
  }

  /**
   * Resets the data and blanks out all charts
   */
  public void clearAllData() {
    ds = new DataStore();
    zooms = ds;
    set = new boolean[DataStore.FILE_COUNT];
    
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
    
    for (ChartPanel cp : chartPanels) {
      JFreeChart chart = ChartFactory.createXYLineChart(
          "",
          "",
          "",
          new XYSeriesCollection(),
          PlotOrientation.VERTICAL,
          false, false, false);
      
      cp.setChart(chart);
      cp.setMouseZoomable(true);
    }
  }
  
}
