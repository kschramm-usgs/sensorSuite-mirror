package asl.sensor;

import java.awt.Dimension;
import java.awt.Insets;
import java.io.File;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;

public class DataPanel extends JPanel {

  /**
   * 
   */
  private static final long serialVersionUID = -7302813951637543526L;

  final static int FILE_COUNT = 3;
  
  // TODO: refactor to be part of InputFileReader?
  TimeSeries[] timeSeriesArray = new TimeSeries[FILE_COUNT];
  
  boolean[] dataSeries = new boolean[FILE_COUNT];

  private JTextArea[] plotters = new JTextArea[FILE_COUNT]; 
              // TODO: replace with chartpanels for data
  
  private ChartPanel[] chartPanels = new ChartPanel[FILE_COUNT];
  
  public DataPanel() {
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    for (int i = 0; i < FILE_COUNT; ++i) {
      
      timeSeriesArray[i] = new TimeSeries("");
      
      JFreeChart chart = ChartFactory.createTimeSeriesChart(
          "",
          "Time",
          "Seismic reading",
          new TimeSeriesCollection(timeSeriesArray[i]),
          false, false, false);
      
      chartPanels[i] = new ChartPanel(chart);
      
      this.add(chartPanels[i]);
      
      /*
      JTextArea statusBox = new JTextArea(5,5);
      statusBox.setMargin( new Insets(5,5,5,5) );
      statusBox.setEditable(false);
      JScrollPane statusScrollPane = new JScrollPane(statusBox);
      statusScrollPane.setBorder( new EmptyBorder(5,5,5,5) );

      
      statusBox.setText("Status Box to eventually hold time-series plots\n");
      statusBox.append("This Box will hold plot for file "+i+'\n');
      plotters[i] = statusBox;
      this.add(statusScrollPane);
      */
      
      this.add( Box.createRigidArea( new Dimension(0,5) ) );
    }
    
  }
  
  public void setData(int idx, String filepath) { 
    // TODO: work with a timeseries object instead of the boolean
    // make sure that chart data is properly updated at that time
    TimeSeries ts = InputFileReader.getTimeSeries(filepath);
    
    timeSeriesArray[idx] = ts;
    
    JFreeChart chart = ChartFactory.createTimeSeriesChart(
        ts.getDomainDescription(),
        "Time",
        "Seismic reading",
        new TimeSeriesCollection(ts),
        false, false, false);
    
    chartPanels[idx].setChart(chart);
    
  }
  
  public TimeSeriesCollection getData() { 
    
    // used to pass data into experimentpanel, to be processed and plotted
    TimeSeriesCollection tsc = new TimeSeriesCollection();
    
    for(TimeSeries timeSeries : timeSeriesArray) {
      tsc.addSeries(timeSeries);
    }
    
    return tsc;
  }
  
  
  
}
