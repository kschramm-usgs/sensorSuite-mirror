package asl.sensor;

public abstract class Experiment {

  // defines template pattern for each type of test, given by backend
  // each test returns new (set of) timeseries data from the input data
  
  // TODO: throw exception if timeSeriesData not same size as FILE_COUNT?
  
  protected double[][] timeSeriesData;
  
  public Experiment(double[][] dataIn) {
    timeSeriesData = backend(dataIn);
  }
  
  public double[][] getData(){
    return timeSeriesData;
  }
  
  abstract double[][] backend(double[][] dataIn);
  
}
