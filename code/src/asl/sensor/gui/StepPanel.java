package asl.sensor.gui;

import java.awt.Color;
import java.awt.Font;

import javax.swing.BoxLayout;

import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

import asl.sensor.experiment.ExperimentEnum;
import asl.sensor.experiment.StepExperiment;
import asl.sensor.input.DataStore;

/**
 * Holds the plot results of a step experiment. Gets the timeseries data from it
 * as well as the corner and damping values gotten in the process of its
 * calculations. See the step experiment class for more details.
 * @author akearns
 *
 */
public class StepPanel extends ExperimentPanel {

  /**
   * 
   */
  private static final long serialVersionUID = 3693391540945130688L;

  /**
   * StepPanel does not have any extra components beyond the abstract panel it
   * is derived from. The constructor initializes its axes and lays out the
   * components derived from the superconstructor
   * @param exp
   */
  public StepPanel(ExperimentEnum exp) {
    super(exp);
    
    xAxisTitle = "Time (s)";
    yAxisTitle = "Normalized counts";
    xAxis = new NumberAxis(xAxisTitle);
    xAxis.setAutoRange(true);
    //SimpleDateFormat sdf = new SimpleDateFormat("HH:mm");
    //sdf.setTimeZone( TimeZone.getTimeZone("UTC") );
    //xAxis.setLabel("UTC Time");
    yAxis = new NumberAxis(yAxisTitle);
    // yAxis.setAutoRange(true);
    Font bold = xAxis.getLabelFont().deriveFont(Font.BOLD);
    xAxis.setLabelFont(bold);
    yAxis.setLabelFont(bold);
    
    applyAxesToChart();
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    this.add(chartPanel);
    this.add(save);
    save.setAlignmentX(CENTER_ALIGNMENT);
    
    plotTheseInBold = new String[]{};
    
    
  }


  /**
   * Pass in and retrieve data from the step experiment backend, to plot;
   * this is both the timeseries data as well as a title inset displaying
   * the parameters used in the plot calculations
   */
  @Override
  public void updateData(final DataStore ds) {
    
    // TODO: threading?
    
    displayInfoMessage("Running stepcal testing...");
    
    expResult.setData(ds, false);
    XYSeriesCollection xysc = expResult.getData();
    
    // here's the stuff that needs to stay here, not moved to experiment class
    populateChart(xysc);
    XYPlot xyp = (XYPlot) chart.getPlot();
    double[] rolloff = ( (StepExperiment) expResult ).getCornerAndDamping();
    double[] fit = ( (StepExperiment) expResult ).getFitCornerAndDamping();
    double corner = rolloff[0];
    double damping = rolloff[1];
    double fitCorner = fit[0];
    double fitDamping = fit[1];
    
    // TODO: will probably need to relocate some of this to its own method
    TextTitle result = new TextTitle();
    StringBuilder sb = new StringBuilder();
    sb.append("RESP parameters\n");
    sb.append("corner frequency: ");
    sb.append(corner);
    sb.append("\n");
    sb.append("damping: ");
    sb.append(damping);
    sb.append("\n");
    sb.append("Best-fit parameters\n");
    sb.append("corner frequency: ");
    sb.append(fitCorner);
    sb.append("\n");
    sb.append("damping: ");
    sb.append(fitDamping);
    sb.append("\n");
    String temp = sb.toString();
    result.setText(temp);
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.5, result,
        RectangleAnchor.RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    chartPanel.setChart(chart);
  }

}