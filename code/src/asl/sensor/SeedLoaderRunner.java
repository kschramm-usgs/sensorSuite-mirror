package asl.sensor;

public class SeedLoaderRunner implements Runnable {

  private String filename;
  private DataBlock blockReadIn;
  
  public SeedLoaderRunner(String filename) {
    this.filename = filename;
  }
  
  public DataBlock getData() {
    return blockReadIn;
  }
  
  @Override
  public void run() {
    // TODO Auto-generated method stub
    blockReadIn = TimeSeriesUtils.getTimeSeries(filename);
  }

}
