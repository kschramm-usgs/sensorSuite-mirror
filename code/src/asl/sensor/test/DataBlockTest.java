package asl.sensor.test;

import static org.junit.Assert.*;

import org.junit.Test;

import asl.sensor.DataBlock;
import asl.sensor.TimeSeriesUtils;
import asl.sensor.DataPanel;

public class DataBlockTest {

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  
  public String fileID = station+"_"+location+"_"+channel;
  
  public String filename1 = "./data/"+fileID+".512.seed";
  
  @Test
  public void trimsCorrectly() {
    int left = 250;
    int right = 750;
    
    DataBlock db = TimeSeriesUtils.getTimeSeries(filename1);

    int sizeOld = db.size();
    
    
    // these get tested in DataPanelTest
    long loc1 = DataPanel.getMarkerLocation(db, left);
    long loc2 = DataPanel.getMarkerLocation(db, right);
    
    db.trim(loc1, loc2);
    
    
    assertEquals( loc1, db.getStartTime() );
    assertEquals( sizeOld/2, db.size() );
  }
  
}