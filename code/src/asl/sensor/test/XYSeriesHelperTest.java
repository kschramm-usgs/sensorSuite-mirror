package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.io.BufferedInputStream;
import java.io.DataInput;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.junit.Test;

import asl.sensor.DataSeriesHelper;
import edu.iris.dmc.seedcodec.B1000Types;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.DecompressedData;
import edu.iris.dmc.seedcodec.UnsupportedCompressionType;
import edu.sc.seis.seisFile.mseed.DataRecord;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import edu.sc.seis.seisFile.mseed.SeedRecord;

public class XYSeriesHelperTest {

  public String station = "TST5";
  public String location = "00";
  public String channel = "BH0";
  
  public String fileID = station+"_"+location+"_"+channel;
  
  public String filename1 = "./data/"+fileID+".512.seed";
  
  @Test
  public void canGetFile() {
    try{
      FileInputStream fis = new FileInputStream(filename1);
      fis.close();
    } catch (Exception e) {
      assertNull(e);
    }
  }
  
  @Test
  public void seisFileCanParseFile() {
    
    try {
      DataInput dis = new DataInputStream( new BufferedInputStream( 
          new FileInputStream(filename1) ) ); 
      try{
        while(true) {
          SeedRecord sr = SeedRecord.read(dis,4096);
          if (sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;

            String loc = dr.getHeader().getLocationIdentifier();
            assertTrue( loc.equals(location) );
            String stat = dr.getHeader().getStationIdentifier().trim();
            assertTrue( stat.equals(station) );

            String chan = dr.getHeader().getChannelIdentifier();
            assertTrue( chan.equals(channel) );
          }
        }
      } catch (EOFException e) {
        assertNotNull(e); // I haaates it! I haaaaaaaaaates it!
      } catch (SeedFormatException e) {
        assertNull(e);
      } catch (IOException e) {
        assertNull(e);
      }
      
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      assertNull(e);
    }
    
  }
  
  @Test
  public void seisFileGivesCorrectSampleRateAndInterval() {
    DataInput dis;
    try {
      while (true) {
        
        dis = new DataInputStream( new BufferedInputStream( 
            new FileInputStream(filename1) ) );
        SeedRecord sr = SeedRecord.read(dis,4096);
        if(sr instanceof DataRecord) {
          DataRecord dr = (DataRecord)sr;
          
          int fact = dr.getHeader().getSampleRateFactor();
          int mult = dr.getHeader().getSampleRateMultiplier();
          
          //System.out.println(fact+","+mult);
          
          double rate = dr.getHeader().getSampleRate();
          assertTrue((double)fact/mult == rate);
          
          // checking the correct values for the intervals
          
          double multOf1Hz = rate/DataSeriesHelper.ONE_HZ;
          long inverse = DataSeriesHelper.ONE_HZ_INTERVAL/(long)multOf1Hz;
          
          long interval = DataSeriesHelper.ONE_HZ_INTERVAL*mult/fact;
          
          assertEquals( inverse, interval);
          // System.out.println(interval);
          
          break;
          
        }
      }
    } catch (FileNotFoundException e) {
      assertNull(e); // only reading one record;
    } catch (SeedFormatException e) {
      assertNull(e);
    } catch (IOException e) {
      assertNull(e);
    }
  }
  
  @Test
  public void inputFileReaderCreatesXYSeries() {
    DataInput dis;
    XYSeries ts = new XYSeries(fileID);
    try {
      dis = new DataInputStream( 
            new BufferedInputStream( 
            new FileInputStream(filename1) ) );
      while ( true ) {
       
        try{
          // going to assume this file's data is continuous. isn't it?
          SeedRecord sr = SeedRecord.read(dis,4096);
          if(sr instanceof DataRecord) {
            
            DataRecord dr = (DataRecord)sr;
            
            int fact = dr.getHeader().getSampleRateFactor();
            int mult = dr.getHeader().getSampleRateMultiplier();
            long interval = DataSeriesHelper.ONE_HZ_INTERVAL*mult/fact;
            
            long startMilli = dr.getHeader()
                                 .getStartBtime()
                                 .convertToCalendar()
                                 .getTimeInMillis();
            
            long activeTime = startMilli;
            
            DecompressedData decomp = dr.decompress();
                        
            // get the original datatype of the series (loads data faster)
            int dataType = decomp.getType();
            // System.out.println(decomp.getTypeString());
            
            // no easy way to do this...
            
            switch (dataType) {
            case B1000Types.INTEGER:
              for (float dataPoint : decomp.getAsInt() ) {
                ts.addOrUpdate(activeTime, dataPoint);
                activeTime += interval;
              }
              break;
            case B1000Types.FLOAT:
              for (float dataPoint : decomp.getAsFloat() ) {
                ts.addOrUpdate(activeTime, dataPoint);
                activeTime += interval;
              }
              break;
            case B1000Types.SHORT:
              for (short dataPoint : decomp.getAsShort() ) {
                ts.addOrUpdate(activeTime, dataPoint);
                activeTime += interval;
              }
              break;
            default: // assume default is double
              for (double dataPoint : decomp.getAsDouble() ) {
                ts.addOrUpdate(activeTime, dataPoint);
                activeTime += interval;
              }
              break;
            }

          }
        
        } catch (EOFException e) {
          break;
        }
      }
      
      XYSeries testAgainst = DataSeriesHelper.getXYSeries(filename1);
      assertEquals( ts.getItemCount(), testAgainst.getItemCount() );
      
    } catch (FileNotFoundException e) {
      assertNull(e);
    } catch (SeedFormatException e) {
      assertNull(e);
    } catch (IOException e) {
      assertNull(e);
    } catch (UnsupportedCompressionType e) {
      // TODO Auto-generated catch block
      assertNull(e);
    } catch (CodecException e) {
      // TODO Auto-generated catch block
      assertNull(e);
    }
  }
  
  
}
