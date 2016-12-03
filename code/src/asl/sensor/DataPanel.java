package asl.sensor;

import java.awt.Dimension;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

public class DataPanel extends JPanel {

  /**
   * 
   */
  private static final long serialVersionUID = -7302813951637543526L;
  
  DataStore ds;
  private ChartPanel[] chartPanels = new ChartPanel[DataStore.FILE_COUNT];
  
  public DataPanel() {
    
    ds = new DataStore();
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    for (int i = 0; i < DataStore.FILE_COUNT; ++i) {
      
      JFreeChart chart = ChartFactory.createTimeSeriesChart(
          "",
          "Time",
          "Seismic reading",
          new TimeSeriesCollection( ds.getSeries(i) ),
          false, false, false);
      
      chartPanels[i] = new ChartPanel(chart);
      //chartPanels[i].setMaximumDrawHeight(50);
      Dimension dim = chartPanels[i].getPreferredSize();
      chartPanels[i].setPreferredSize(
          new Dimension( (int) dim.getWidth(), (int) dim.getHeight()/2));
      chartPanels[i].setMouseZoomable(false);
      
      this.add(chartPanels[i]);
      
      this.add( Box.createRigidArea( new Dimension(0,5) ) );
    }
    
  }
  
  public void setData(int idx, String filepath) { 
    
    TimeSeries ts = ds.setData(idx, filepath);
    
    JFreeChart chart = ChartFactory.createTimeSeriesChart(
        ts.getKey().toString(),
        "Time",
        "Seismic reading",
        new TimeSeriesCollection(ts),
        false, false, false);
    
    chartPanels[idx].setChart(chart);
    chartPanels[idx].setMouseZoomable(false);
    
  }
  
  public TimeSeriesCollection getData() { 
    
    return ds.getData();

  }
  
  
  
}
