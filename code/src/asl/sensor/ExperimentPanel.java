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

import org.apache.commons.math3.complex.Complex;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeries;
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
  public static JFreeChart populateChart(ExperimentEnum expT, Experiment expR,
      boolean freqSpace) {
    
    XYSeriesCollection data = expR.getData();
    
    String xTitle;
    if (freqSpace) {
      xTitle = expR.getXTitle();
    } else {
      xTitle = expR.getFreqTitle();
    }
    
    JFreeChart chart = ChartFactory.createXYLineChart(
        expT.getName(),
        xTitle,
        expR.getYTitle(),
        data,
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    
    XYPlot xyPlot = chart.getXYPlot();
    
    String[] bold = expR.getBoldSeriesNames();
    XYItemRenderer xyir = xyPlot.getRenderer();
    for (String series : bold) {
      int seriesIdx = data.getSeriesIndex(series);
      BasicStroke stroke = (BasicStroke) xyir.getBaseStroke();
      stroke = new BasicStroke( stroke.getLineWidth()*2 );
      xyir.setSeriesStroke(seriesIdx, stroke);
      xyir.setSeriesPaint(seriesIdx, new Color(0,0,0) );
    }
    
    if (freqSpace) {
      xyPlot.setDomainAxis(expR.getFreqAxis());
    } else {
      xyPlot.setDomainAxis(expR.getXAxis());

    }
    xyPlot.setRangeAxis(expR.getYAxis());
    
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
    chartPanel.setMouseZoomable(false);
    
    fc = new JFileChooser();
    
    save = new JButton("Save Plot (PNG)");
    save.addActionListener(this);
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
  
  /**
   * Used to propagate changes in data into the underlying experiment process.
   * @param tsc The new data, in DataBlock format 
   *            (List, startTime, name, interval [inverse of sample rate])
   */
  public void updateData(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    expResult.setData(ds, psd, freqSpace);
    
    chart = populateChart(expType, expResult, freqSpace);
    
    chartPanel.setChart(chart);
    chartPanel.setMouseZoomable(false);
    
    // setting the new chart is enough to update the plots
  }
  
  public static void addToPlot(
      DataBlock[] dataIn, 
      FFTResult[] psd, 
      XYSeriesCollection plottable,
      boolean freqSpace) {
    
    for (int i = 0; i < psd.length; ++i) {

      XYSeries powerSeries = new XYSeries( "PSD " + dataIn[i].getName() );

      Complex[] resultPSD = psd[i].getFFT();
      double[] freqs = psd[i].getFreqs();

      for (int j = 0; j < freqs.length; ++j) {
        if (1/freqs[j] > 1.0E3) {
          continue;
        }

        // TODO: is this right (seems to be)
        Complex temp = resultPSD[j].multiply(Math.pow(2*Math.PI*freqs[j],4));

        if (freqSpace) {
          powerSeries.add( freqs[j], 10*Math.log10( temp.abs() ) );
        } else {
          powerSeries.add( 1/freqs[j], 10*Math.log10( temp.abs() ) );
        }
      }

      plottable.addSeries(powerSeries);

    }
    
  }

  public abstract void setDataNames(String[] seedFileNames);
  
}
