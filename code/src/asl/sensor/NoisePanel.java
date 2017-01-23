package asl.sensor;

import java.awt.Checkbox;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.text.SimpleDateFormat;
import java.util.TimeZone;

import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.SwingWorker;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class NoisePanel extends ExperimentPanel {

  private Checkbox freqSpaceBox;
  
  public NoisePanel(ExperimentEnum exp) {
    // create chart, chartPanel, save button & file chooser, 
    super(exp);
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    // chart = populateChart(expType, expResult);
    
    freqSpaceBox = new Checkbox("Use Hz units (requires regen)");
    freqSpaceBox.setState(false);
    
    JPanel optionsPanel = new JPanel();
    optionsPanel.setLayout( new BoxLayout(optionsPanel, BoxLayout.Y_AXIS) );
    optionsPanel.add(freqSpaceBox);
    optionsPanel.add(save);
    save.setAlignmentX(CENTER_ALIGNMENT);

    
    // do the layout here
    this.add(chartPanel);
    this.add(optionsPanel);

    
  }
  
  @Override
  public void actionPerformed(ActionEvent e) {
    super.actionPerformed(e); // only actionlistener here
  }
  
  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    
    if (ds.numberFullySet() < 3) {
      displayErrorMessage("INSUFFICIENT DATA LOADED");
      return;
    }
    
    boolean freqSpace = freqSpaceBox.getState();
    
    updateDriver(ds, psd, freqSpace);
    // setting the new chart is enough to update the plots
    
    
  }
  
  private void updateDriver(DataStore ds, FFTResult[] psd, boolean freqSpace) {
    
    final DataStore dsImmutable = ds;
    final FFTResult[] psdImmutable = psd;
    final boolean freqSpaceImmutable = freqSpace;
    
    SwingWorker<Integer, Void> worker = new SwingWorker<Integer, Void>() {
      @Override
      public Integer doInBackground() {
        expResult.setData(dsImmutable, psdImmutable, freqSpaceImmutable);
        
        chart = populateChart(expResult.getData());
        
        // override the default axis if the checkbox is set to use Hz units
        if (freqSpaceImmutable) {
          ValueAxis xAxis = 
              ( (NoiseExperiment) expResult ).getXAxis(freqSpaceImmutable);
          chart.getXYPlot().setDomainAxis(xAxis);
        }
        return 0;
      }

      @Override
      public void done() {
        chartPanel.setChart(chart);
        chartPanel.setMouseZoomable(false);
      }

    };
    
    new Thread(worker).start();
    return;
    
  }
  

  /**
   * Auto-generated serialize ID
   */
  private static final long serialVersionUID = 9018553361096758354L;
  
  

}
