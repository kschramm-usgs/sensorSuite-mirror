package asl.sensor.test;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.junit.Test;

import asl.sensor.InstrumentResponse;
import asl.sensor.TransferFunction;
import asl.sensor.Unit;

public class InstrumentResponseTest {

  @Test
  public void testFileParse() {
    String filename = "./responses/RESP.XX.NS087..BHZ.STS1.20.2400";
    
    try {
      InstrumentResponse ir = new InstrumentResponse(filename);
      
      assertEquals(TransferFunction.LAPLACIAN, ir.getTransferFunction() );
      
      double nml = Double.parseDouble("3.948580E+03");
      assertEquals( nml, ir.getNormalization(), 0.0001 );
      
      double nmf = Double.parseDouble("3.000000E-01");
      assertEquals( nmf, ir.getNormalizationFrequency(), 0.0001 );
      
      double[] gn = {2.400000e+03, 2.400000e+03, 1.000000e+00};
      assertTrue( Arrays.equals( gn, ir.getGain() ) );
      
      assertEquals( Unit.VELOCITY, ir.getUnits() );
      
      List<Complex> zrs = new ArrayList<Complex>();
      zrs.add( new Complex(0.000000e+00, 0.000000e+00) );
      zrs.add( new Complex(0.000000e+00, 0.000000e+00) );
      assertEquals( zrs, ir.getZeros() );
      
      List<Complex> pls = new ArrayList<Complex>();
      pls.add( new Complex(-2.221000e-01,  2.221000e-01) );
      pls.add( new Complex(-2.221000e-01, -2.221000e-01) );
      pls.add( new Complex(-3.918000e+01,  4.912000e+01) );
      pls.add( new Complex(-3.918000e+01, -4.912000e+01) );
      assertEquals( pls, ir.getPoles() );
      
    } catch (IOException e) {
      // TODO Auto-generated catch block
      fail("Unexpected error trying to read response file");
    }
    
  }

}
