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

public class TimeSeriesHelper {
  
  // TODO: 2 slow 4 me, need to fix that

  public final static long ONE_HZ_INTERVAL = 1000000L;
  public final static double ONE_HZ = 1.0;
  
  public final static int MAX_SIZE = 1000000;

  public static TimeSeries reduce(TimeSeries ts) {
    if( ts != null &&  ts.getItemCount() > TimeSeriesHelper.MAX_SIZE) {
      int plot_size = ts.getItemCount()/2;
      int factor = 2;
      while (plot_size > TimeSeriesHelper.MAX_SIZE) {
        plot_size /= 2;
        factor *= 2;
      }
      
      TimeSeries shorter = new TimeSeries(ts.getKey().toString());
      
      for (int i = 0; i < ts.getItemCount(); i+=factor ){
        shorter.add( ts.getDataItem(i) );
      }

      return shorter;
    } else {
      return ts;
    }
    
  }
  
  
  public static TimeSeries getTimeSeries(String filename) {
    
    // TODO: easiest way to speed this up would be to
    // filter out points not within a range of interest BEFORE plotting, etc.
    // Does any such method of extracting data exist?
    
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