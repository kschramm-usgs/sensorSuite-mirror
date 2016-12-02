package asl.sensor;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.RegularTimePeriod;
import org.jfree.data.time.TimeSeries;

import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.UnsupportedCompressionType;
import edu.sc.seis.seisFile.mseed.DataHeader;
import edu.sc.seis.seisFile.mseed.DataRecord;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import edu.sc.seis.seisFile.mseed.SeedRecord;

public class InputFileReader {

  public static final long ONE_HZ_INTERVAL = 1000000L;
  public static final double ONE_HZ = 1.0;
  
  // the seisfile library has a convenient DataRecord object type
  // that can make it easy to get SEED file components like sample rate,
  // station name, etc. and of course the sweet juicy data candy in the center
  
  // important notes:
  
  // basic structure can probably be taken from miniseedtofloatarray example
  
  // the time interval can be taken from sampleratefactor and sampleratemult.
  // both of which are publicly accessible from a datarecord's getters
  
  // add the first sample to a timeseries object at the starting time long
  // then the next sample's time is that plus the interval (1/sampleRate)
  // TODO: see what sample rate factor and multiple actually are
  // (want to avoid as many conversions to/from float due to noise)
  
  public static TimeSeries getTimeSeries(String filename) {
    
    DataInputStream dis;
    TimeSeries ts = null;
    
    try {
      dis = new DataInputStream(  new FileInputStream(filename) );

      while ( true ) {
        
        try {
          SeedRecord sr = SeedRecord.read(dis,4096);
          if(sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;
            DataHeader dh = dr.getHeader();
            if (ts == null){
              StringBuilder fileID = new StringBuilder();
              fileID.append(dh.getStationIdentifier() + "_");
              fileID.append(dh.getLocationIdentifier() + "_");
              fileID.append(dh.getChannelIdentifier());
              ts = new TimeSeries(fileID.toString());
            }
            
            // TODO: do we need a check that file datarecords are all contiguous?
            long start = dh.getStartBtime()
                           .convertToCalendar()
                           .getTimeInMillis();
            long active = start;
            
            int fact = dh.getSampleRateFactor();
            int mult = dh.getSampleRateMultiplier();
            
            long interval = ONE_HZ_INTERVAL*mult/fact;
            
            double[] data = dr.decompress().getAsDouble();
            
            // TODO: can this be made faster?
            
            for(double dataPoint : data) {
              RegularTimePeriod sampleTaken = new FixedMillisecond(active); 
              ts.addOrUpdate(sampleTaken, dataPoint);
              active += interval;
            }
          }
        } catch(EOFException e) {
          break;
        }
          
      }
      
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (SeedFormatException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (UnsupportedCompressionType e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (CodecException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    return ts;
  }
  

}
