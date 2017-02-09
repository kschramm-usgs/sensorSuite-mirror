package asl.sensor;

import java.awt.Color;
import java.awt.Font;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;

import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class StepPanel extends ExperimentPanel {

  /**
   * 
   */
  private static final long serialVersionUID = 3693391540945130688L;
  
  JFileChooser fc;

  public StepPanel(ExperimentEnum exp) {
    super(exp);
    
    xAxisTitle = "Time (s)";
    yAxisTitle = "Counts (normalized)";
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

  @Override
  public void updateData(DataStore ds) {
    
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
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.02, result,
        RectangleAnchor.BOTTOM_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    chartPanel.setChart(chart);
  }

}
