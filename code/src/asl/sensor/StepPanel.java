package asl.sensor;

import java.util.Arrays;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JPanel;

import org.apache.commons.math3.complex.Complex;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.JLabel;

public class StepPanel extends ExperimentPanel {

  JFileChooser fc;
  JButton loadCal;
  JLabel fileNameArea;
  
  
  public StepPanel(ExperimentEnum exp) {
    super(exp);
    
    fc = new JFileChooser();
    
    this.setLayout( new BoxLayout(this, BoxLayout.Y_AXIS) );
    
    loadCal = new JButton("Load calibration data");
    fileNameArea = new JLabel("NO FILE LOADED");
    
    this.add(chartPanel);
    
    JPanel loadPanel = new JPanel();
    loadPanel.setLayout( new BoxLayout(loadPanel, BoxLayout.X_AXIS) );
    loadPanel.add(loadCal);
    loadPanel.add(fileNameArea);
    this.add(loadPanel);
    
  }

  /**
   * 
   */
  private static final long serialVersionUID = 3693391540945130688L;

  @Override
  public void updateData(DataStore ds, FFTResult[] psd) {
    // TODO Auto-generated method stub
    
    // TODO move to stepexperiment
    DataBlock db = ds.getBlock(0);
    long interval = db.getInterval();
    InstrumentResponse ir = ds.getResponse(0);
    Complex pole = ir.getPoles().get(0);
    double f = 1. / (2 * Math.PI / pole.abs() ); // corner frequency
    double h = Math.abs( pole.getReal() / pole.abs() ); // damping
    
    System.out.println(f+","+h);
    
    InstrumentResponse appxResp = new InstrumentResponse(f, h);
    
    // get FFT of datablock timeseries, apply response to input
    FFTResult fft = FFTResult.simpleFFT(db);
    
    Complex[] respValues = appxResp.applyResponseToInput( fft.getFreqs() );
    
    double maxVal = respValues[0].abs();
    for (Complex respVal : respValues) {
      if ( respVal.abs()  > maxVal ) {
        maxVal = respVal.abs();
      }
    }
    
    System.out.println(maxVal);
    
    Complex[] fftValues = fft.getFFT();
    
    Complex[] correctedValues = new Complex[fftValues.length];
    // don't let denominator be zero
    for (int i = 0; i < correctedValues.length; ++i) {
      Complex numer = fftValues[i].multiply( respValues[i].conjugate() );
      Complex denom = respValues[i].multiply( respValues[i].conjugate() );
      denom = denom.add(0.008 * maxVal);
      correctedValues[i] = numer.divide(denom);
    }
    
    double[] toPlot = FFTResult.inverseFFT(correctedValues);
    long now = 0;
    
    XYSeries xys = new XYSeries( db.getName() );
    for (double point : toPlot) {
      xys.add(point, now);
      now += interval;
    }
    
    
    XYSeriesCollection xysc = new XYSeriesCollection();
    xysc.addSeries(xys);
    
    // here's the stuff that needs to stay here, not moved to experiment class
    chart = populateChart(xysc, true);
    chartPanel.setChart(chart);
  }

}
