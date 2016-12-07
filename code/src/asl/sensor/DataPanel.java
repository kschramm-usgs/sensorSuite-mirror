package asl.sensor;

import java.awt.Color;
import java.awt.Component;
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
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


/**
 * Panel used to hold the plots for the files taken in as input
 * @author akearns
 *
 */
public class DataPanel extends JPanel implements ActionListener {

  /**
   * 
   */
  private static final long serialVersionUID = -7302813951637543526L;
  
  DataStore ds;
  private ChartPanel[] chartPanels = new ChartPanel[DataStore.FILE_COUNT];
  private Color[] defaultColor = {
          ChartColor.LIGHT_RED, 
          ChartColor.LIGHT_BLUE, 
          ChartColor.LIGHT_GREEN };
  private JButton save;
  private JFileChooser fc;
  private JPanel allCharts; // parent of the chartpanels, used for image saving
  
  
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
          "Time",
          "Seismic reading",
          new XYSeriesCollection( ds.getPlotSeries(i) ),
          PlotOrientation.VERTICAL,
          false, false, false);
      
      chartPanels[i] = new ChartPanel(chart);
      //chartPanels[i].setMaximumDrawHeight(50);
      Dimension dim = chartPanels[i].getPreferredSize();
      chartPanels[i].setPreferredSize(
          new Dimension( (int) dim.getWidth(), (int) dim.getHeight()/2));
      chartPanels[i].setMouseZoomable(true);
      
      allCharts.add(chartPanels[i]);
      
      // don't add a space below the last plot (yet)
      if( i+1 < DataStore.FILE_COUNT) {
        allCharts.add( Box.createRigidArea( new Dimension(0,5) ) );
      }

      
    }
    
    this.add(allCharts);
    
    // now we can add the space between the last plot and the save button
    this.add( Box.createRigidArea( new Dimension(0,5) ) );
    
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
        "Seismic reading",
        new XYSeriesCollection(ts),
        PlotOrientation.VERTICAL,
        false, false, false);
    
    XYPlot xyp = (XYPlot) chart.getPlot();
    xyp.getRenderer().setSeriesPaint(0, defaultColor[idx]);
    
    chartPanels[idx].setChart(chart);
    chartPanels[idx].setMouseZoomable(true);
    
  }
  
  /**
   * Returns the underlying DataStore's image, to be fed into experiments
   * for processing (the results of which will be plotted)
   * @return An array of DataBlocks (time series and metadata), indexed by plot
   *         order.
   */
  public DataBlock[] getData() { 
    
    return ds.getData();

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
          
          BufferedImage bi = new BufferedImage(
                  allCharts.getWidth(), 
                  allCharts.getHeight(), 
                  BufferedImage.TYPE_INT_ARGB);
          
          Graphics2D g = bi.createGraphics();
          allCharts.printAll(g);
          g.dispose();
          
          ImageIO.write(bi,"png",selFile);
        } catch (IOException e1) {
          // TODO Auto-generated catch block
          e1.printStackTrace();
        }
      }
    }
  }
  
  
  
}
