package asl.sensor;

import org.jfree.chart.annotations.XYTitleAnnotation;
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
    xAxisTitle = "Ortho X Axis (units)";
    yAxisTitle = "Ortho Y Axis (units)";
    // TODO Auto-generated constructor stub
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
    sb.append(fit[0]);
    sb.append('\n');
    sb.append(fit[1]);
    TextTitle result = new TextTitle();
    result.setText( sb.toString() );
    XYTitleAnnotation xyt = new XYTitleAnnotation(0.98, 0.98, result,
        RectangleAnchor.TOP_RIGHT);
    xyp.clearAnnotations();
    xyp.addAnnotation(xyt);
    
    chartPanel.setChart(chart);
    
  }

}
