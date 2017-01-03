package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
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
   * 
   */
  private static final long serialVersionUID = -7302813951637543526L;
  
  public static final int IMAGE_HEIGHT = 240*DataStore.FILE_COUNT;
  
  DataStore ds;
  private ChartPanel[] chartPanels = new ChartPanel[DataStore.FILE_COUNT];
  private Color[] defaultColor = {
          ChartColor.LIGHT_RED, 
          ChartColor.LIGHT_BLUE, 
          ChartColor.LIGHT_GREEN };
  private JButton save;
  private JFileChooser fc;
  private JPanel allCharts; // parent of the chartpanels, used for image saving
  private JPanel sliderPanel;
  private JSlider leftSlider;
  private JSlider rightSlider;
  private int rightSliderValue;
  private int leftSliderValue;
  private int margin = 100; // min space of the two sliders
  
  
  /**
   * Creates a new data panel -- instantiates each chart, to be populated with
   * data when a file is loaded in. Also creates a save button for writing all
   * the inputted data plots into a single PNG file.
   */
  public DataPanel() {
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
   
    ds = new DataStore();
    
    allCharts = new JPanel();
    allCharts.setLayout( new BoxLayout(allCharts, BoxLayout.Y_AXIS) );
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      
      JFreeChart chart = ChartFactory.createXYLineChart(
          "",
          "",
          "",
          new XYSeriesCollection( ds.getPlotSeries(i) ),
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
    
    leftSliderValue = 0;
    rightSliderValue = 1000;
    
    sliderPanel = new JPanel();
    sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.X_AXIS));
    sliderPanel.setBorder(new EmptyBorder(5, 10, 5, 10));
    leftSlider = new JSlider(0, 1000, leftSliderValue);
    leftSlider.addChangeListener(this);
    sliderPanel.add(leftSlider);
    rightSlider = new JSlider(0, 1000, rightSliderValue);
    rightSlider.addChangeListener(this);
    sliderPanel.add(rightSlider);
    this.add(sliderPanel);
    
    // now we can add the space between the last plot and the save button
    this.add( Box.createVerticalStrut(5) );
    
    save = new JButton("Save");
    this.add(save);
    save.setAlignmentX(Component.CENTER_ALIGNMENT);
    save.addActionListener(this);
    
    fc = new JFileChooser();
    
  }
  
  /**
   * Takes a loaded and converted miniSEED time series and plots it
   * by calling the underlying DataStore
   * @param idx Index of chart to be loaded to (0 to DataStore.FILE_COUNT)
   * @param filepath The full address of the file to be loaded in
   */
  public void setData(int idx, String filepath) { 
    
    XYSeries ts = ds.setData(idx, filepath);
    
    JFreeChart chart = ChartFactory.createXYLineChart(
        ts.getKey().toString(),
        "Time",
        "Counts",
        new XYSeriesCollection(ts),
        PlotOrientation.VERTICAL,
        false, false, false);
    
    XYPlot xyp = (XYPlot) chart.getPlot();
    DateAxis da = new DateAxis();
    SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");
    sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    da.setLabel("UTC Time");
    Font bold = da.getLabelFont();
    bold = bold.deriveFont(Font.BOLD);
    da.setLabelFont(bold);
    da.setDateFormatOverride(sdf);
    xyp.setDomainAxis(da);
    xyp.getRenderer().setSeriesPaint(0, defaultColor[idx]);
    
    chartPanels[idx].setChart(chart);
    chartPanels[idx].setMouseZoomable(true);
    
    leftSlider.setValue(0);
    rightSlider.setValue(1000);
    setVerticalBars();
    
  }
  
  /**
   * Displays the range set by the sliders using
   * vertical bars at the min and max values
   */
  public void setVerticalBars() {
    
    ds.trimToCommonTime();
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      if ( null == ds.getBlock(i) ) {
        leftSlider.setValue(0);
        rightSlider.setValue(1000);
        return;
      }
      
      int leftValue = leftSlider.getValue();
      int rightValue = rightSlider.getValue();
      
      XYPlot xyp = (XYPlot) chartPanels[i].getChart().getPlot();
      xyp.clearDomainMarkers();
      
      DataBlock db = ds.getBlock(i);
      
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
    
  }
  
  /**
   * Checks if all the data has been loaded
   * @return True if all data is set, including responses
   */
  public boolean dataIsSet() {
    return ds.isPlottable();
  }
  
  
  /**
   * Returns the selected region of underlying DataStore, to be fed into experiments
   * for processing (the results of which will be plotted)
   * @return A DataStore object (contains arrays of DataBlocks & Responses)
   */
  public DataStore getData() {
    if ( ! ds.isPlottable() ){
      throw new RuntimeException("Not all necessary data loaded in...");
    }
    // TODO: modify this to create a new dataStore with trimmed data
    int leftValue = leftSlider.getValue();
    int rightValue = rightSlider.getValue();
    DataBlock db = ds.getBlock(0); // dataBlocks should have same time range
    long start = getMarkerLocation(db, leftValue);
    long end = getMarkerLocation(db, rightValue);
    return new DataStore(ds, start, end);
  }

  /**
   * Handles saving this panel's plots to file (PNG image)
   * when the save button is clicked.
   */
  @Override
  public void actionPerformed(ActionEvent e) {

    if( e.getSource() == save ) {
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
    }
  }
  
  /**
   * Return this panel's charts as a single buffered image
   * @return Buffered image of the plots, writeable to file
   */
  public BufferedImage getAsImage() {
    return getAsImage( allCharts.getWidth(), allCharts.getHeight() );
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
    
    width = Math.max( width, chartPanels[0].getWidth() );
    // height = Math.max( height, shownHeight );
    
    // cheap way to make sure height is a multiple of the chart count
    height = (height*DataStore.FILE_COUNT)/DataStore.FILE_COUNT;
    
    int chartHeight = height/DataStore.FILE_COUNT;
    
    Dimension outSize = new Dimension(width, height);

    JPanel toDraw = new JPanel(); // what we're going to draw to image
    toDraw.setSize(outSize);
    toDraw.setPreferredSize(outSize);
    toDraw.setMinimumSize(outSize);
    toDraw.setMaximumSize(outSize);
    
    toDraw.setLayout( new BoxLayout(toDraw, BoxLayout.Y_AXIS) );

    for (ChartPanel cp : chartPanels) {
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
  
  
}
