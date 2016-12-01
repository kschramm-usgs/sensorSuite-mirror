package asl.sensor;

import java.awt.Component;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingConstants;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.general.SeriesException;
import org.jfree.data.time.Millisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

public class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;
  // likely need a JFreeChart object of some kind here

  private JButton save;
  private JFreeChart chart; // replace with plot object
  
  private TimeSeriesCollection datasets;
  
  // private JTextArea txa;
  
  private Experiment exp_; // used to define experiment of each plot object
    // TODO: replace with an experiment plotter once available
  
  public ExperimentPanel(Experiment exp){
    
    exp_ = exp;
    
    createDataset();
    
    chart = ChartFactory.createXYLineChart(
        exp_.getName(),
        "X",
        "Y",
        datasets);
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    save = new JButton("Save Plot");
    save.addActionListener(this);

    /*
    txa = new JTextArea(20,50);
    txa.setMargin( new Insets(5,5,5,5) );
    txa.setEditable(false);
    txa.setText("GUI PLOT AREA WIP "+exp.getName()+'\n');
    
    // just to make testing easier
    txa.append( txa.getSize().getHeight() + "" );
    */
    

    
    //this.add(textScroll);
    this.add( new ChartPanel(chart) );
    this.add(save);
    save.setAlignmentX(Component.CENTER_ALIGNMENT);
    
  }
  
  private void createDataset(){
    datasets = new TimeSeriesCollection();
    
    Millisecond step = new Millisecond();
    double value = 100.0;
    for(int i = 0; i < DataPanel.FILE_COUNT; ++i) {
      for (int j = 0; j < 4000; ++j) {
         TimeSeries series = new TimeSeries("Dataset " + i);
         try {
            value = value + Math.random( ) - 0.5;                 
            series.add( step, new Double(value) );                 
            step = ( Millisecond ) step.next( ); 
         } catch ( SeriesException e ) {
            System.err.println("Error adding to series");
         }
      }
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
