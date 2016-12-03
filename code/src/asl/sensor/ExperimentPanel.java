package asl.sensor;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.time.TimeSeriesCollection;

public class ExperimentPanel extends JPanel implements ActionListener {
 
  private static final long serialVersionUID = -5591522915365766604L;
  // likely need a JFreeChart object of some kind here

  private JButton save;
  private JFreeChart chart; // replace with plot object
  private ChartPanel chartPanel;
  
  private JFileChooser fc; // save image when image save button clicked
  
//  private TimeSeriesCollection datasets;
  
  private ExperimentEnum expType; 
          // used to define experiment of each plot object
  
  private Experiment expResult;
          // used to get the actual data from loaded-in files
  
  public static JFreeChart populateChart(ExperimentEnum expT, Experiment expR) {
    JFreeChart chart = ChartFactory.createTimeSeriesChart(
        expT.getName(),
        expR.getXTitle(),
        expR.getYTitle(),
        expR.getData(),
        true, // include legend
        false, 
        false);
    
    return chart;
  }
  
  public ExperimentPanel(ExperimentEnum exp, TimeSeriesCollection tsc) {
    
    expType = exp;
    // TODO: reset with actual input data
    expResult = ExperimentFactory.createExperiment(exp, tsc);
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    save = new JButton("Save Plot");
    save.addActionListener(this);

    chartPanel = new ChartPanel( populateChart(expType, expResult) );
    chartPanel.setMouseZoomable(false);
    
    this.add(chartPanel);
    this.add(save);
    save.setAlignmentX(Component.CENTER_ALIGNMENT);
    
    fc = new JFileChooser();
    
  }
  
  public void updateData(TimeSeriesCollection tsc) {
    
    expResult.setData(tsc);
    
    chartPanel.setChart( populateChart(expType, expResult) );
    chartPanel.setMouseZoomable(false);
    
    // setting the new chart is enough to update the plots
  }

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
          ChartUtilities.saveChartAsPNG(selFile,chart,1280,960);
        } catch (IOException e1) {
          // TODO Auto-generated catch block
          e1.printStackTrace();
        }
      }
    }
  }
  
}
