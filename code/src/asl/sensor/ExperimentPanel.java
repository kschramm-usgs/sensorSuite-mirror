package asl.sensor;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.general.SeriesException;
import org.jfree.data.time.Second;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

public class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;
  // likely need a JFreeChart object of some kind here

  private JButton save;
  private JFreeChart chart; // replace with plot object
  
  private TimeSeriesCollection datasets;
  
  private Experiment exp_; // used to define experiment of each plot object
    // TODO: replace with an experiment plotter once available
  
  public ExperimentPanel(Experiment exp){
    
    exp_ = exp;
    
    createDataset();
    
    chart = ChartFactory.createTimeSeriesChart(
        exp_.getName(),
        "Time (sec)",
        "Value",
        datasets,
        true, // include legend
        false, 
        false);
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    save = new JButton("Save Plot");
    save.addActionListener(this);

    ChartPanel chartp = new ChartPanel(chart);
    chartp.setMouseZoomable(false);
    
    this.add(chartp);
    this.add(save);
    save.setAlignmentX(Component.CENTER_ALIGNMENT);
    
  }
  
  private void createDataset(){
    datasets = new TimeSeriesCollection();
    
    for(int i = 0; i < DataPanel.FILE_COUNT; ++i) {
      Second step = new Second(0,0,0,1,1,2016);
      double value = 100.0;
      TimeSeries series = new TimeSeries("Dataset " + i);
      for (int j = 0; j < 4000; ++j) {
         try {
            value = value + Math.random( ) - 0.5;                 
            series.add( step, new Double(value) );                 
            step = ( Second ) step.next( );
            System.out.println("flow test, " + value + " " + step);
         } catch ( SeriesException e ) {
            System.err.println("Error adding to series");
         }
      }
      datasets.addSeries(series);
    }
  }

  @Override
  public void actionPerformed(ActionEvent e) {
    // TODO Auto-generated method stub
    if( e.getSource() == save ) {
      // TODO: save plot to file
      // txa.append("\nClicked the save button.");
    }
    
  }
  
}
