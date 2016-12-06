package asl.sensor;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.jfree.data.xy.XYSeries;

import edu.iris.dmc.seedcodec.B1000Types;
import edu.iris.dmc.seedcodec.CodecException;
import edu.iris.dmc.seedcodec.DecompressedData;
import edu.iris.dmc.seedcodec.UnsupportedCompressionType;
import edu.sc.seis.seisFile.mseed.DataHeader;
import edu.sc.seis.seisFile.mseed.DataRecord;
import edu.sc.seis.seisFile.mseed.SeedFormatException;
import edu.sc.seis.seisFile.mseed.SeedRecord;

public class DataSeriesHelper {
  
  // TODO: 2 slow 4 me, need to fix that

  public final static long ONE_HZ_INTERVAL = 1000000L;
  public final static double ONE_HZ = 1.0;
  
  public static XYSeries getXYSeries(String filename) {
    
    // TODO: what can we do to make this faster?
    
    DataInputStream dis;
    XYSeries xys = null;
    
    try {
      dis = new DataInputStream(  new FileInputStream(filename) );

      while ( true ) {
        
        try {
          SeedRecord sr = SeedRecord.read(dis,4096);
          if(sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;
            DataHeader dh = dr.getHeader();
            if (xys == null){
              StringBuilder fileID = new StringBuilder();
              fileID.append(dh.getStationIdentifier() + "_");
              fileID.append(dh.getLocationIdentifier() + "_");
              fileID.append(dh.getChannelIdentifier());
              xys = new XYSeries(fileID.toString());
            }
            
            
            byte af = dh.getActivityFlags();
            byte correctionFlag = 0b00000010; // is there a time correction?
            byte qf = dh.getDataQualityFlags();
            byte glitchFlag = 0b00001000; // does the record have errors?
            if( ( qf & glitchFlag) != 0) {
              continue;
            }
            int correction = 0;
            if ( (af & correctionFlag) == 0 ) {
              correction = dh.getTimeCorrection();
            }

            long start = dh.getStartBtime()
                           .convertToCalendar()
                           .getTimeInMillis() + correction;
            long active = start;
            
            int fact = dh.getSampleRateFactor();
            int mult = dh.getSampleRateMultiplier();
            
            long interval; // mult inverse of sample rate
            
            if( fact > 0 && mult > 0) {
              interval = ONE_HZ_INTERVAL / (fact * mult);
            } else if (fact > 0 && mult < 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * mult) / fact);
            } else if (fact < 0 && mult > 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * fact) / mult);
            } else {
              interval = ONE_HZ_INTERVAL * fact * mult;
            }
            
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
                xys.addOrUpdate(active, dataPoint);
                active += interval;
              }
              break;
            case B1000Types.FLOAT:
              for (float dataPoint : decomp.getAsFloat() ) {
                xys.addOrUpdate(active, dataPoint);
                active += interval;
              }
              break;
            case B1000Types.SHORT:
              for (short dataPoint : decomp.getAsShort() ) {
                xys.addOrUpdate(active, dataPoint);
                active += interval;
              }
              break;
            default:
              for (double dataPoint : decomp.getAsDouble() ) {
                xys.addOrUpdate(active, dataPoint);
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
    
    return xys;
  }
  

}
