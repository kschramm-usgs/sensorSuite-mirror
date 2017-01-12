package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Panel used to display the data produced from a specified sensor test
 * @author akearns
 *
 */
public abstract class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;

  /**
   * Defines the type of test results this panel's chart will display
   * @param expT The type of experiment (refer to the enum for valid types)
   * @param expR The experiment to be performed
   * @return A chart displaying the data from the performed experiment
   */
  public JFreeChart populateChart(XYSeriesCollection data, boolean freqSpace) {
    
    String xTitle;
    if (freqSpace) {
      xTitle = expResult.getXTitle();
    } else {
      xTitle = expResult.getFreqTitle();
    }
    
    JFreeChart chart = ChartFactory.createXYLineChart(
        expType.getName(),
        xTitle,
        expResult.getYTitle(),
        data,
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    
    XYPlot xyPlot = chart.getXYPlot();
    
    String[] bold = expResult.getBoldSeriesNames();
    XYItemRenderer xyir = xyPlot.getRenderer();
    for (String series : bold) {
      int seriesIdx = data.getSeriesIndex(series);
      BasicStroke stroke = (BasicStroke) xyir.getBaseStroke();
      stroke = new BasicStroke( stroke.getLineWidth()*2 );
      xyir.setSeriesStroke(seriesIdx, stroke);
      xyir.setSeriesPaint(seriesIdx, new Color(0,0,0) );
    }
    
    if (freqSpace) {
      xyPlot.setDomainAxis(expResult.getFreqAxis());
    } else {
      xyPlot.setDomainAxis(expResult.getXAxis());

    }
    xyPlot.setRangeAxis(expResult.getYAxis());
    
    return chart;
  }
  
  protected JButton save;
  protected JFreeChart chart; // replace with plot object
  
  protected ChartPanel chartPanel;
  
  protected JFileChooser fc; // save image when image save button clicked
  
  protected ExperimentEnum expType; 
          // used to define experiment of each plot object
  
  protected Experiment expResult;
          // used to get the actual data from loaded-in files
  
  public ExperimentPanel(ExperimentEnum exp) {
    
    expType = exp;
    expResult = ExperimentFactory.createExperiment(exp);
    
    chart = ChartFactory.createXYLineChart(expType.getName(), 
        expResult.getXTitle(), expResult.getYTitle(), null);
    chartPanel = new ChartPanel(chart);
    // chartPanel.setMouseZoomable(false);
    
    fc = new JFileChooser();
    
    save = new JButton("Save Plot (PNG)");
    save.addActionListener(this);
  }
  
  public boolean isTwoInput() {
    return expType.isTwoInput();
  }
  
  /**
   * Handle's saving this plot's chart to file (PNG image) 
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
          ChartUtilities.saveChartAsPNG(selFile,chart,640,480);
        } catch (IOException e1) {
          e1.printStackTrace();
        }
      }
    }
  }

  /**
   * Return image of this chart with specified dimensions
   * Used to compile PNG image of all currently-displayed charts
   * @param width Width of output image in pixels
   * @param height Height of output image in pixels
   * @return buffered image of this panel's chart
   */
  public BufferedImage getAsImage(int width, int height) {
    
    ChartPanel outPanel = new ChartPanel(chart);
    outPanel.setSize( new Dimension(width, height) );
    
    BufferedImage bi = new BufferedImage(
        width, 
        height, 
        BufferedImage.TYPE_INT_ARGB);

    Graphics2D g = bi.createGraphics();
    outPanel.printAll(g);
    g.dispose();

    return bi;
  }
  
  public abstract void updateData(DataStore ds, FFTResult[] psd);

  public void updateData(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    expResult.setData(ds, psd, freqSpace);
    
    chart = populateChart(expResult.getData(), freqSpace);
    
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);
    
    // setting the new chart is enough to update the plots
  }

  public void setDataNames(String[] seedFileNames) {
    // nothing to use the filenames for here
    return;
  }
  
}
