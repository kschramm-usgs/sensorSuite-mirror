package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.TimeZone;

import javax.imageio.ImageIO;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.SwingWorker;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;

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
 * @author akearns
 *
 */
public class DataPanel 
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
  private JFileChooser fc;
  private JPanel allCharts; // parent of the chartpanels, used for image saving
  private JPanel sliderPanel;
  private JSlider leftSlider;
  private JSlider rightSlider;
  private int margin = 100; // min space of the two sliders
  
  
  /**
   * Creates a new data panel -- instantiates each chart, to be populated with
   * data when a file is loaded in. Also creates a save button for writing all
   * the inputted data plots into a single PNG file.
   */
  public DataPanel() {
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
   
    ds = new DataStore();
    zooms = ds;
    
    allCharts = new JPanel();
    allCharts.setLayout( new BoxLayout(allCharts, BoxLayout.Y_AXIS) );
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      
      JFreeChart chart = ChartFactory.createXYLineChart(
          "",
          "",
          "",
          new XYSeriesCollection(),
          PlotOrientation.VERTICAL,
          false, false, false);
      
      chartPanels[i] = new ChartPanel(chart);
      Dimension dim = chartPanels[i].getPreferredSize();
      chartPanels[i].setPreferredSize(
          new Dimension( (int) dim.getWidth(), (int) dim.getHeight()/2) );
      chartPanels[i].setMaximumSize(
          new Dimension(Integer.MAX_VALUE, Integer.MAX_VALUE) );
      chartPanels[i].setMouseZoomable(true);
      
      allCharts.add(chartPanels[i]);
      
      // don't add a space below the last plot (yet)
      if( i+1 < DataStore.FILE_COUNT) {
        allCharts.add( Box.createVerticalStrut(5) );
      }

      
    }
    
    this.add(allCharts);
    
    sliderPanel = new JPanel();
    sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.X_AXIS));
    sliderPanel.setBorder(new EmptyBorder(5, 10, 5, 10));
    leftSlider = new JSlider(0, 1000, 0);
    leftSlider.setEnabled(false);
    leftSlider.addChangeListener(this);
    sliderPanel.add(leftSlider);
    rightSlider = new JSlider(0, 1000, 1000);
    rightSlider.setEnabled(false);
    rightSlider.addChangeListener(this);
    sliderPanel.add(rightSlider);
    this.add(sliderPanel);
    
    // now we can add the space between the last plot and the save button
    this.add( Box.createVerticalStrut(5) );
    
    JPanel buttons = new JPanel();
    buttons.setLayout( new BoxLayout(buttons, BoxLayout.X_AXIS) );
    
    zoomIn = new JButton("Zoom on Selection");
    buttons.add(zoomIn);
    zoomIn.addActionListener(this);
    zoomIn.setEnabled(false);
    
    save = new JButton("Save Timeseries Plots (PNG)");
    buttons.add(save);
    save.addActionListener(this);
    
    zoomOut = new JButton("Show all data");
    buttons.add(zoomOut);
    zoomOut.addActionListener(this);
    zoomOut.setEnabled(false);
    
    this.add(buttons);
    
    fc = new JFileChooser();
    
  }
  
  
  /**
   * Takes a loaded and converted miniSEED time series and plots it
   * by calling the underlying DataStore
   * @param idx Index of chart to be loaded to (0 to DataStore.FILE_COUNT)
   * @param filepath The full address of the file to be loaded in
   */
  public void setData(int idx, String filepath) { 

    // set all data to the same range first (zoom out)
    zooms = ds;
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      if (zooms.getBlock(i) == null) {
        continue;
      }
      this.resetPlotZoom(i);
    }
    
    final String pathFinal = filepath;
    final int index = idx;

    SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {

      JFreeChart chart;
      
      @Override
      public Integer doInBackground() {
        ds.setData(index, pathFinal);
        XYSeries ts = ds.getBlock(index).toXYSeries();
        try {
          chart = ChartFactory.createXYLineChart(
              ts.getKey().toString(),
              "Time",
              "Counts",
              new XYSeriesCollection(ts),
              PlotOrientation.VERTICAL,
              false, false, false);

          XYPlot xyp = (XYPlot) chart.getPlot();
          DateAxis da = new DateAxis();
          SimpleDateFormat sdf = new SimpleDateFormat("HH:mm");
          sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
          da.setLabel("UTC Time");
          Font bold = da.getLabelFont();
          bold = bold.deriveFont(Font.BOLD);
          da.setLabelFont(bold);
          da.setDateFormatOverride(sdf);
          xyp.setDomainAxis(da);
          xyp.getRenderer().setSeriesPaint(0, defaultColor[index]);

          return 0;
        } catch (Exception e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        return 1;
      }
      
      @Override
      public void done() {
        // drawing the chart generally takes longer than reading in the data
        // which can be very fast if we don't use the error checking read-in
        // so threading here may not be necessary
        
        if (chart == null) {
          set[index] = false;
          return;
        }
        
        chartPanels[index].setChart(chart);
        chartPanels[index].setMouseZoomable(true);

        showRegionForGeneration();
        
        set[index] = true;

        zoomIn.setEnabled(true);
        leftSlider.setEnabled(true);
        rightSlider.setEnabled(true);

        leftSlider.setValue(0);
        rightSlider.setValue(1000);
        setVerticalBars();
      }
      
    };

    new Thread(worker).start();
    return;

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
    
    int leftValue = leftSlider.getValue();
    int rightValue = rightSlider.getValue();
    DataBlock db = zooms.getBlock(0); // default initialization
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      // dataBlocks should have same time range
      db = zooms.getBlock(i); 
      if (db == null) {
        // if this data wasn't set, try the next one
        continue;
      } else {
        break; // all data should have the same range
      }
    }

    long start = getMarkerLocation(db, leftValue);
    long end = getMarkerLocation(db, rightValue);
    return new DataStore(zooms, start, end);
  }

  /**
   * Dispatches commands based on save and zoom buttons clicked
   */
  @Override
  public void actionPerformed(ActionEvent e) {

    if ( e.getSource() == save ) {
      String ext = ".png";
      fc.addChoosableFileFilter(
          new FileNameExtensionFilter("PNG image (.png)",ext) );
      fc.setFileFilter(fc.getChoosableFileFilters()[1]);
      int returnVal = fc.showSaveDialog(save);
      if (returnVal == JFileChooser.APPROVE_OPTION) {
        File selFile = fc.getSelectedFile();
        if( !selFile.getName().endsWith( ext.toLowerCase() ) ) {
          selFile = new File( selFile.toString() + ext);
        }
        try {
          // each chart gets displayed as 640x240 image; merge to one file
          BufferedImage bi = getAsImage(640, IMAGE_HEIGHT);
          
          ImageIO.write(bi,"png",selFile);
        } catch (IOException e1) {
          // TODO Auto-generated catch block
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
    int height = allCharts.getHeight();
    // otherwise, we only use the height of images actually set
    height /= DataStore.FILE_COUNT;
    height *= ds.numberOfBlocksLoaded();
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
    // TODO: don't bother including plots that have no data loaded
    // (replace FILE_COUNT with a call to 'amountOfDataLoaded' or similar)
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
          leftSliderValue + margin > rightSliderValue) {
        leftSliderValue = rightSliderValue - margin;
        if (leftSliderValue < 0) {
          leftSliderValue = 0;
          rightSliderValue = margin;
        }
      }
    } else if ( e.getSource() == rightSlider ) {
      if (rightSliderValue < leftSliderValue ||
          rightSliderValue - margin < leftSliderValue) {
        rightSliderValue = leftSliderValue + margin;
        if (rightSliderValue > 1000) {
          rightSliderValue = 1000;
          leftSliderValue = 1000-margin;
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
