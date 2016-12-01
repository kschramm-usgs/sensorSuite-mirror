package asl.sensor;

public abstract class Experiment {

  double[][] timeSeriesData;  
  
  abstract void backend(double[][] dataIn);
  
}
