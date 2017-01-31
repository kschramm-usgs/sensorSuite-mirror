package asl.sensor.test;

import static org.junit.Assert.*;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import org.junit.Test;

import asl.sensor.DataBlock;
import asl.sensor.TimeSeriesUtils;
import asl.sensor.DataPanel;

public class DataPanelTest {

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  
  public String fileID = station+"_"+location+"_"+channel;
  
  public String filename1 = "./data/"+fileID+".512.seed";
  
  @Test
  public void getsCorrectTrimming() {
    int left = 250;
    int right = 750;
    
    String name;
    try {
      name = new ArrayList<String>( 
          TimeSeriesUtils.getMplexNameSet(filename1)
        ).get(0);
      
      DataBlock db = TimeSeriesUtils.getTimeSeries(filename1, name);
      long start = db.getStartTime();
      long interval = db.getInterval();
      int size = db.size();
      
      long timeRange = interval*size;
      
      long loc1 = DataPanel.getMarkerLocation(db, left);
      long loc2 = DataPanel.getMarkerLocation(db, right);
      
      assertEquals(loc2 - loc1, timeRange/2); // 500/1000 = 1/2
      assertEquals(loc1, start + (interval * 1/4 * size) ); // correct start pt?
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
      fail();
    }
  
    
    
  }
  
}
