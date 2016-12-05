package asl.sensor;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.jfree.data.time.FixedMillisecond;
import org.jfree.data.time.RegularTimePeriod;
import org.jfree.data.time.TimeSeries;

import edu.iris.dmc.seedcodec.B1000Types;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.DecompressedData;
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
    
    // TODO: what can we do to make this faster?
    
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
            
            DecompressedData decomp = dr.decompress();

            // get the original datatype of the series (loads data faster)
            // otherwise the decompressed data gets converted (cloned) as
            // the other type instead
            int dataType = decomp.getType();
            
            // This is probably the best way to do this since
            // we have to add each point individually anyway
            // and converting between types for 

            switch (dataType) {
            case B1000Types.INTEGER:
              for (int dataPoint : decomp.getAsInt() ) {
                RegularTimePeriod thisTimeStamp = new FixedMillisecond(active);
                ts.addOrUpdate(thisTimeStamp, dataPoint);
                active += interval;
              }
              break;
            case B1000Types.FLOAT:
              for (float dataPoint : decomp.getAsFloat() ) {
                RegularTimePeriod thisTimeStamp = new FixedMillisecond(active);
                ts.addOrUpdate(thisTimeStamp, dataPoint);
                active += interval;
              }
              break;
            case B1000Types.SHORT:
              for (short dataPoint : decomp.getAsShort() ) {
                RegularTimePeriod thisTimeStamp = new FixedMillisecond(active);
                ts.addOrUpdate(thisTimeStamp, dataPoint);
                active += interval;
              }
              break;
            default:
              for (double dataPoint : decomp.getAsDouble() ) {
                RegularTimePeriod thisTimeStamp = new FixedMillisecond(active);
                ts.addOrUpdate(thisTimeStamp, dataPoint);
                active += interval;
              }
              break;
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
