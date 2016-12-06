package asl.sensor;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

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
  
  public static DataBlock getXYSeries(String filename) {
    
    // TODO: what can we do to make this faster?
    
    DataInputStream dis;
    // XYSeries xys = null;
    DataBlock db = null;
    List<Number> data = new ArrayList<Number>();
    
    try {
      dis = new DataInputStream(  new FileInputStream(filename) );

      while ( true ) {
        
        try {
          long interval = 0L;
          SeedRecord sr = SeedRecord.read(dis,4096);
          if(sr instanceof DataRecord) {
            DataRecord dr = (DataRecord)sr;
            DataHeader dh = dr.getHeader();
            if (db == null){
              StringBuilder fileID = new StringBuilder();
              fileID.append(dh.getStationIdentifier() + "_");
              fileID.append(dh.getLocationIdentifier() + "_");
              fileID.append(dh.getChannelIdentifier());
              db = new DataBlock(data, interval, fileID.toString(), -1);
            }
            
            
            byte af = dh.getActivityFlags();
            byte correctionFlag = 0b00000010; // is there a time correction?
            int correction = 0;
            if ( (af & correctionFlag) == 0 ) {
              correction = dh.getTimeCorrection();
            }

            long start = dh.getStartBtime()
                           .convertToCalendar()
                           .getTimeInMillis() + correction;
            
            if(db.startTime < 0) {
              db.startTime = start;
            }
            
            int fact = dh.getSampleRateFactor();
            int mult = dh.getSampleRateMultiplier();
               
            if( fact > 0 && mult > 0) {
              interval = ONE_HZ_INTERVAL / (fact * mult);
            } else if (fact > 0 && mult < 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * mult) / fact);
            } else if (fact < 0 && mult > 0) {
              interval = Math.abs( (ONE_HZ_INTERVAL * fact) / mult);
            } else {
              interval = ONE_HZ_INTERVAL * fact * mult;
            }
            
            db.interval = interval;
            
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
              int[] decomArrayInt = decomp.getAsInt();
              for (int dataPoint : decomArrayInt ) {
                data.add(dataPoint);
              }
              break;
            case B1000Types.FLOAT:
              float[] decomArrayFlt = decomp.getAsFloat();
              for (float dataPoint : decomArrayFlt ) {
                data.add(dataPoint);
              }
              break;
            case B1000Types.SHORT:
              short[] decomArrayShr = decomp.getAsShort();
              for (short dataPoint : decomArrayShr ) {
                data.add(dataPoint);
              }
              break;
            default:
              double[] decomArrayDbl = decomp.getAsDouble();
              for (double dataPoint : decomArrayDbl ) {
                data.add(dataPoint);
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
    
    db.data = data;
    
    return db;
  }
  
  /**
   * Initial driver for the decimation utility
   * which takes a timeseries of unknown rate and
   * runs downsampling to convert it to a target
   * frequency of a 1Hz interval.
   * @param data The timeseries to be decimated
   * @param src The source frequency as interval between samples (microseconds)
   * @return A timeseries decimated to the correct frequency
   */
  public static List<Number> decimate(List<Number> data, long src){
     
     long tgt = ONE_HZ_INTERVAL; // target frequency
     // a sample lower than 1Hz frq has longer time between samples
     // since it's an inverse relationship and all
     if(src >= tgt){
       // if data is too low-frequency to decimate, do nothing
       return data;
     }

     // find what the change in size is going to be
     long gcd = euclidGCD(src, tgt);
     // conversion up- and down-factors
     // (upsample by target, downsample by source)
     // cast is valid because any discrete interval
     // from 1Hz and up is already expressable
     // as an int
     int upf = (int)(src/gcd);
     int dnf = (int)(tgt/gcd);

     // one valid sample rate for data is 2.5Hz
     // with 1Hz that comes out as a ratio of 5/2, which won't
     // downsample neatly in some cases so we would first upsample,
     // filter out any noise terms, then downsample
     List<Number> upped = upsample(data,upf);
     List<Number> lpfed = lowPassFilter(upped,src*upf);
     List<Number> down = downsample(lpfed,dnf);
     
     return down;

  }

  /**
   * Implements Euclid's algorithm for finding GCD
   * used to find common divisors to give us upsample
   * and downsample rates by dividing the timeseries intervals
   * by this value
   * @param src Initially, one of two frequencies to calculate
   * @param tgt Initially, one of two frequencies to calculate
   * @return The GCD of the two frequencies
   */
  public static long euclidGCD(long src,long tgt){
  
    // take remainders until we hit 0
    // which means the divisor is the gcd
    long rem = src % tgt;
    if(rem == 0){
      return tgt;
    }
    
    return euclidGCD(tgt, rem);
  }

  /**
   * Upsamples data by a multiple of passed factor, placing zeros
   * between each data point. Result is data.length*factor cells in size.
   * Requires use of a low-pass filter to remove discontinuities.
   * @param data The timeseries to be upsampled
   * @param factor The factor to increase the size by
   * @return The upsampled series
   */
  public static List<Number> upsample(List<Number> data, int factor){

    List<Number> upsamp = Arrays.asList(new Number[data.size()*factor]);
   
    for(int i=0; i<data.size(); i++){
      upsamp.set( i*factor, data.get(i) ); // index, element
    }
    
    return upsamp;
  }

  /**
   * Downsamples data by a multiple of passed factor. Result is
   * data.length/factor cells in size
   * Requires previous use of a low-pass filter to avoid aliasing
   * @param data The timeseries to be downsampled
   * @param factor The factor to decrease the size by
   * @return The downsampled series
   */
  public static List<Number> downsample(List<Number> data, int factor){
    
    List<Number> downsamp = Arrays.asList(new Number[data.size()/factor]);
    for(int i=0; i < downsamp.size(); i++){
      downsamp.set( i, data.get(i*factor) ); 
    }
    
    return downsamp;
  }

  /**
   * Implements low pass band filter
   * @param timeseries  The data to be filtered
   * @param sps         Samples per second
   * @return            The filtered data
   */
  public static List<Number> lowPassFilter(List<Number> timeseries, long sps)
  {
    // apache fft requires input to be power of two
    int pow = 2;
    while( pow < timeseries.size() ){
      pow *= 2;
    }
    
    double[] timeseriesFilter = new double[pow];

    List<Number> timeseriesdouble = 
                  Arrays.asList( new Number[timeseries.size()] );
    double fl = 5; // allow all low-frequency data through
    double fh = ONE_HZ_INTERVAL/2.0; // nyquist rate half of target frequency
       // note that this 1/2F where F is 1Hz frequency
    // we want the inverse of the target frequency

    for (int ind = 0; ind < timeseries.size(); ind++)
    {
      if( !(timeseries.get(ind) instanceof Double) ) {
        timeseriesFilter[ind] = (double) timeseries.get(ind).intValue();
      }
    }

    FastFourierTransformer fft = 
                  new FastFourierTransformer(DftNormalization.STANDARD);
    
    
    Complex[] frqDomn = fft.transform(timeseriesFilter, TransformType.FORWARD);
    
    frqDomn = apply((double) sps, frqDomn, fl, fh);

    Complex[] timeDomn = fft.transform(frqDomn, TransformType.INVERSE);
    

    for (int ind = 0; ind < timeseries.size(); ind++)
    {
      timeseriesdouble.set(ind, (double) timeDomn[ind].getReal() );
    }

    return timeseriesdouble;
  }

  /**
   * Implements bandpass filter for lowpassfilter()
   * @param dt  Time step
   * @param cx  Complex number form of time series
   * @param fl  Low corner frequency
   * @param fh  High corner frequency
   * @return  Complex form of filtered time series
   */
  public static Complex[] apply(double dt, Complex[] cx, double fl, double fh)
  {

    int npts = cx.length;
    // double fl = 0.01;
    // double fh = 2.0;
    int npole = 2;
    int numPoles = npole;
    int twopass = 2;
    double TWOPI = Math.PI * 2;
    double PI = Math.PI;

    Complex c0 = new Complex(0., 0.);
    Complex c1 = new Complex(1., 0.);

    Complex[] sph = new Complex[numPoles];
    Complex[] spl = new Complex[numPoles];

    Complex cjw, cph, cpl;
    int nop, nepp, np;
    double wch, wcl, ak, ai, ar, w, dw;
    int i, j;

    if (npole % 2 != 0)
    {
      System.out.println("WARNING - Number of poles not a multiple of 2!");
    }

    nop = npole - 2 * (npole / 2);
    nepp = npole / 2;
    wch = TWOPI * fh;
    wcl = TWOPI * fl;

    np = -1;
    if (nop > 0)
    {
      np = np + 1;
      sph[np] = new Complex(1., 0.);
    }
    if (nepp > 0)
    {
      for (i = 0; i < nepp; i++)
      {
        ak = 2. * Math
           .sin((2. * (double) i + 1.0) * PI / (2. * (double) npole));
        ar = ak * wch / 2.;
        ai = wch * Math.sqrt(4. - ak * ak) / 2.;
        np = np + 1;
        sph[np] = new Complex(-ar, -ai);
        np = np + 1;
        sph[np] = new Complex(-ar, ai);
      }
    }
    np = -1;
    if (nop > 0)
    {
      np = np + 1;
      spl[np] = new Complex(1., 0.);
    }
    if (nepp > 0)
    {
      for (i = 0; i < nepp; i++)
      {
        ak = 2. * Math
           .sin((2. * (double) i + 1.0) * PI / (2. * (double) npole));
        ar = ak * wcl / 2.;
        ai = wcl * Math.sqrt(4. - ak * ak) / 2.;
        np = np + 1;
        spl[np] = new Complex(-ar, -ai);
        np = np + 1;
        spl[np] = new Complex(-ar, ai);
      }
    }

    cx[0] = c0;
    dw = TWOPI / ((double) npts * dt);
    w = 0.;
    for (i = 1; i < npts / 2 + 1; i++)
    {
      w = w + dw;
      cjw = new Complex(0., -w);
      cph = c1;
      cpl = c1;
      for (j = 0; j < npole; j++)
      {
        
        cph = cph.multiply(sph[j]);
        Complex div = new Complex(0.0, 0.0);
        div = sph[j].add(cjw);
        cph.divide(div);
        
        cpl = cpl.multiply(cjw);
        div = new Complex(0.0, 0.0);
        div = spl[j].add(cjw);
        cpl = cpl.divide(div);
        
        /*
        cph = Complex.div(Complex.mul(cph, sph[j]), Complex.add(sph[j], cjw));
        cpl = Complex.div(Complex.mul(cpl, cjw), Complex.add(spl[j], cjw));
        */
      }
      Complex prod = new Complex(0.0,0.0).add(cph).multiply(cpl).conjugate();
      cx[i] = cx[i].multiply(prod);

      // cx[i] = Complex.mul(cx[i], (Complex.mul(cph, cpl)).conjg());
      
      if (twopass == 2)
      {
        cx[i] = cx[i].multiply( cph.add(cpl) );
      }
      cx[npts - i] = (cx[i]).conjugate();
    }

    return (cx);

  }

}

