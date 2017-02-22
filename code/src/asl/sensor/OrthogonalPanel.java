package asl.sensor;

import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class OrthogonalPanel extends ExperimentPanel {

  /**
   * 
   */
  private static final long serialVersionUID = -2749224338484110043L;

  public OrthogonalPanel(ExperimentEnum exp) {
    super(exp);
    xAxisTitle = "Time (s)";
    yAxisTitle = "Amplitude difference (counts)";
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
    
    this.setLayout( new GridBagLayout() );
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = 0; gbc.gridy = 0;
    gbc.weightx = 1.0; gbc.weighty = 1.0;
    gbc.fill = GridBagConstraints.BOTH;
    gbc.anchor = GridBagConstraints.CENTER;
    this.add(chartPanel, gbc);
    gbc.weighty = 0.0;
    gbc.fill = GridBagConstraints.NONE;
    gbc.gridy += 1;
    this.add(save, gbc);
    
  }

  @Override
  public void updateData(DataStore ds) {
    // TODO Auto-generated method stub
    
    expResult.setData(ds, false);
    XYSeriesCollection xysc = expResult.getData();
    
    populateChart(xysc);
    XYPlot xyp = (XYPlot) chart.getPlot();
    
    double[] fit = ( (OrthogonalExperiment) expResult ).getSolutionParams();
    
    StringBuilder sb = new StringBuilder();
    sb.append("Parameters for fit relative to reference sensor:\n");
    sb.append(fit[0] + " (rel. reference #1)");
    sb.append('\n');
    sb.append(fit[1] + " (rel. reference #2)");
    TextTitle result = new TextTitle();
    result.setText( sb.toString() );
    result.setBackgroundPaint(Color.white);
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    chartPanel.setChart(chart);
    
  }

}
