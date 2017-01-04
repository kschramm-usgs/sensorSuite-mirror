package asl.sensor.test;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import asl.sensor.DataBlock;
import asl.sensor.DataPanel;
import asl.sensor.DataStore;

public class DataStoreTest {

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  
  public String fileID = station+"_"+location+"_"+channel;
  
  public String filename1 = "./data/"+fileID+".512.seed";
  
  @Test
  public void commonTimeTrimMatchesLength() {
    
    DataStore ds = new DataStore();
    ds.setData(0,filename1);
    
    int left = 250;
    int right = 750;
    
    DataBlock db = ds.getBlock(0);
    
    // tested in DataPanelTest
    long loc1 = DataPanel.getMarkerLocation(db, left);
    long loc2 = DataPanel.getMarkerLocation(db, right);
    //  tested in DataBlockTest
    db.trim(loc1, loc2);
    
    ds.setData(1,filename1);
    ds.setData(2,filename1);
    
    // function under test
    ds.trimToCommonTime();
    
    assertEquals( ds.getBlock(1).getStartTime(), loc1);
    assertEquals( db.size(), ds.getBlock(1).size() );
    
  }
  
}