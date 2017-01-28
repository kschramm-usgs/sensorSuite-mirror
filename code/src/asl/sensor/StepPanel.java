package asl.sensor;

import java.awt.Color;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;

import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleAnchor;

public class StepPanel extends ExperimentPanel {

  JFileChooser fc;
  
  public StepPanel(ExperimentEnum exp) {
    super(exp);
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    this.add(chartPanel);
    this.add(save);
    save.setAlignmentX(CENTER_ALIGNMENT);
    
  }

  /**
   * 
   */
  private static final long serialVersionUID = 3693391540945130688L;

  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    // TODO Auto-generated method stub
    
    int numLoaded = ds.numberOfBlocksLoaded();
    int numLoadedWithResp = ds.numberFullySet();
    
    // we need two sets of data, one of which has a response file
    if (numLoaded < 2 || numLoadedWithResp < 1) {
      displayErrorMessage("INSUFFICIENT DATA LOADED");
      return;
    }

    
    expResult.setData(ds, psd, false);
    XYSeriesCollection xysc = expResult.getData();
    
    // here's the stuff that needs to stay here, not moved to experiment class
    chart = populateChart(xysc);
    XYPlot xyp = (XYPlot) chart.getPlot();
    double[] rolloff = ( (StepExperiment) expResult ).getCornerAndDamping();
    double corner = rolloff[0];
    double damping = rolloff[1];
    
    // TODO: will probably need to relocate some of this to its own method
    TextTitle result = new TextTitle();
    StringBuilder sb = new StringBuilder();
    sb.append("corner frequency: ");
    sb.append(corner);
    sb.append("\n");
    sb.append("damping: ");
    sb.append(damping);
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
