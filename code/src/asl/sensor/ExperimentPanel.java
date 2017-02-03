package asl.sensor;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GridBagLayout;
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
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

/**
 * Panel used to display the data produced from a specified sensor test
 * @author akearns
 *
 */
public abstract class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;

  protected JButton save;
  
  protected JFreeChart chart; // replace with plot object
  protected ChartPanel chartPanel;
  
  protected JFileChooser fc; // save image when image save button clicked
  
  protected ExperimentEnum expType; 
          // used to define experiment of each plot object
  
  protected Experiment expResult;
          // used to get the actual data from loaded-in files
  
  protected String xAxisTitle, yAxisTitle;
  protected ValueAxis xAxis, yAxis;
  protected String[] plotTheseInBold;
  
  public ExperimentPanel(ExperimentEnum exp) {
    
    expType = exp;
    expResult = ExperimentFactory.createExperiment(exp);
    
    chart = ChartFactory.createXYLineChart(expType.getName(), 
        "", "", null);
    chartPanel = new ChartPanel(chart);
    // chartPanel.setMouseZoomable(false);
    
    fc = new JFileChooser();
    
    save = new JButton("Save plot (PNG)");
    save.addActionListener(this);
    
    this.setLayout( new GridBagLayout() );
  }
  
  /**
   * Since the constructor here must call the 
   */
  protected void applyAxesToChart() {
    XYPlot xyp = chart.getXYPlot();
    xyp.setDomainAxis( getXAxis() );
    xyp.setRangeAxis( getYAxis() );
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
  
  public void displayInfoMessage(String infoMsg) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText(infoMsg);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
  }
  
  public void displayErrorMessage(String errMsg) {
    XYPlot xyp = (XYPlot) chartPanel.getChart().getPlot();
    TextTitle result = new TextTitle();
    result.setText(errMsg);
    result.setBackgroundPaint(Color.red);
    result.setPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.5, 0.5, result,
        RectangleAnchor.CENTER);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
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

  public boolean isTwoInput() {
    return expType.isTwoInput();
  }
  

  public String getXTitle() {
    return xAxisTitle;
  }

  public String getYTitle() {
    return yAxisTitle;
  }


  public ValueAxis getXAxis() {
    // TODO Auto-generated method stub
    return xAxis;
  }


  public ValueAxis getYAxis() {
    // TODO Auto-generated method stub
    return yAxis;
  }

  
  public String[] seriesToDrawBold() {
    return plotTheseInBold;
  }
  
  
  /**
   * Defines the type of test results this panel's chart will display
   * @param expT The type of experiment (refer to the enum for valid types)
   * @param expR The experiment to be performed
   * @return A chart displaying the data from the performed experiment
   */
  public void populateChart(XYSeriesCollection data) {
    
    String xTitle = getXTitle();
    
    JFreeChart chart = ChartFactory.createXYLineChart(
        expType.getName(),
        xTitle,
        getYTitle(),
        data,
        PlotOrientation.VERTICAL,
        true, // include legend
        false, 
        false);
    
    // apply bolding to the components that require it (i.e., NLNM time series)
    XYPlot xyPlot = chart.getXYPlot();
    String[] bold = seriesToDrawBold();
    XYItemRenderer xyir = xyPlot.getRenderer();
    for (String series : bold) {
      int seriesIdx = data.getSeriesIndex(series);
      BasicStroke stroke = (BasicStroke) xyir.getBaseStroke();
      stroke = new BasicStroke( stroke.getLineWidth()*2 );
      xyir.setSeriesStroke(seriesIdx, stroke);
      xyir.setSeriesPaint(seriesIdx, new Color(0,0,0) );
    }
    
    xyPlot.setDomainAxis( getXAxis() );
    xyPlot.setRangeAxis( getYAxis() );
    
    this.chart = chart;
  }
  
  public abstract void updateData(DataStore ds);
  
  // details of how to run updateData are left up to the implementing panel
  // however, it is advised to wrap the code inside a swingworker,
  // using doInBackground() to make calls to the experiment backend to
  // compute the data and then using the done() method to update chartPanel

}
