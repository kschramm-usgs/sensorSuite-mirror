package asl.sensor;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.complex.Complex;

public class InstrumentResponse {

  public static final int GAIN_STAGES = 3; // sensitivity, stage 1, stage 2
  
  private TransferFunction transferType;
  
  // gain values, indexed by stage
  private double[] gain = new double[GAIN_STAGES];
  
  private List<Complex> zeros;
  private List<Complex> poles;
  
  private Unit unitType;
  
  private double normalization;
  private double normalFreq; // ♫ cuz she's a normalFreq, normalFreq 
  
  public InstrumentResponse(String filename) throws IOException {
    parseResponseFile(filename);
  }
  
  public InstrumentResponse(InstrumentResponse responseIn) {
    transferType = responseIn.getTransferFunction();
    
    int len = responseIn.getGain().length;
    gain = Arrays.copyOf( responseIn.getGain(), len );

    zeros = new ArrayList<Complex>( responseIn.getZeros() );
    poles = new ArrayList<Complex>( responseIn.getPoles() );
    
    unitType = responseIn.getUnits();
    
    normalization = Double.valueOf( responseIn.getNormalization() );
    normalFreq = Double.valueOf( responseIn.getNormalizationFrequency() );
  }
  
  public TransferFunction getTransferFunction() {
    return transferType;
  }
  
  public double[] getGain() {
    return gain;
  }
  
  public List<Complex> getZeros() {
    return zeros;
  }
  
  public List<Complex> getPoles() {
    return poles;
  }
  
  public Unit getUnits() {
    return unitType;
  }
  
  public double getNormalization() {
    return normalization;
  }
  
  public double getNormalizationFrequency() {
    return normalFreq;
  }
  
  public Complex[] applyResponseToInput(double[] frequencies) {
   
    Complex[] resps = new Complex[frequencies.length];
    
    // precalculate gain for scaling the response
    // use gain1 * gain2 unless gain0 differs by more than 10%
    // (apparently an issue with Q680 detectors)
    double diff = 100 * ( gain[0] - (gain[1]*gain[2]) ) / gain[0];
    double scale;
    if (Math.abs(diff) > 10) {
      scale = gain[0];
    } else {
      scale = gain[1] * gain[2];
    }
    
    
    for (int i =0; i < frequencies.length; ++i) {
      // TODO: check this is right (Adam seemed to be saying so)
      double deltaFrq = frequencies[i];
      
      // pole-zero expansion
      Complex s = new Complex( 0, deltaFrq*transferType.getFunction() );
      
      Complex numerator = Complex.ONE;
      Complex denominator = Complex.ONE;
      
      for (Complex zero : zeros) {
        numerator = numerator.multiply( s.subtract(zero) );
      }
      
      for (Complex pole : poles) {
        denominator = denominator.multiply( s.subtract(pole) );
      }
      
      resps[i] = numerator.multiply(normalization).divide(denominator);
      
      // how many times do we need to do differentiation?
      int differentiations = Unit.ACCELERATION.getDifferentiations(unitType);
      
      // TODO: also do a check for integration?
      for (int j = 0; j < differentiations; ++j ) {
        // unlike S this is always 2Pi
        double integConstant = 2*Math.PI;
        // i*omega -- differentiate n times by taking I(omg) / (i*omg)^n
        Complex iw = new Complex(0.0, -1.0 / (integConstant * frequencies[i]) );
        // (just need to do this one to go from vel to accel)
        resps[i] = iw.multiply(resps[i]);
      }
      
      // lastly, scale by the scale we chose (gain0 or gain1*gain2)
      resps[i] = resps[i].multiply(scale);
      
    }
    
    return resps;
  }
  
  /**
   * Parses a response file of the sort found on the Iris Nominal Response
   * Library. These files can be found at http://ds.iris.edu/NRL/
   * This function currently does not parse a full response file, but instead
   * only examines fields relevant to self-noise calculations.
   * @param filename Full path to the response file
   */
  private void parseResponseFile(String filename) throws IOException {
    
    // response files have a very nice format that is not so nice as something
    // like JSON but still quite easy to parse
    // lines either begin with a hex value or a '#'
    // '#' marks comments while the hex value is used to interpret values
    // each line with a hex value appears to have the following format
    // [hex value] [whitespace] [name] [whitespace] [value]
    // where name is the human-readable explanation of what a value represents
    // in some cases, 'value' may not be just a raw value but also include
    // some information about the value, such as verbose unit specifications
    
    // there is one exception, the actual pole/zero fields, which have 5
    // components after the hex identifier
    
    BufferedReader br;
    try {
      br = new BufferedReader( new FileReader(filename) );
      String line = br.readLine();
      
      int gainStage = -1;
      Complex[] polesArr = null;
      Complex[] zerosArr = null;
      
      while (line != null) {
        if (line.charAt(0) == '#') {
          // comment -- skip
          line = br.readLine();
          continue;
        } else {
          // the components of each line, assuming split by 2 or more spaces
          String[] words = line.split("[ ]{2,}");
          String hexIdentifier = words[0];
          
          switch (hexIdentifier) {
          case "B053F03":
            // transfer function type specified
            // first character of third component of words
            switch ( words[2].charAt(0) ) {
            case 'A':
              transferType = TransferFunction.LAPLACIAN;
              break;
            case 'B':
              transferType = TransferFunction.LINEAR;
              break;
            default:
              // TODO: want to make sure that we're setting this right
              // because if a Fourier is specified, we should set to LAPLACIAN
              transferType = TransferFunction.LAPLACIAN;
            }
            break;
          case "B053F05":
            // parse the units of the transfer function (usually velocity)
            // first *word* of the third component of words
            // TODO support more unit types
            String[] unitString = words[2].split(" ");
            String unit = unitString[0];
            switch (unit.toLowerCase()) {
            case "m/s":
              unitType = Unit.VELOCITY;
              break;
            case "m/s**2":
              unitType = Unit.ACCELERATION;
              break;
            default:
              String e = "Nonstandard units, or not a velocity or acceleration";
              throw new IOException(e);
            }
            break;
          case "B053F07":
            // this is the normalization factor A0
            // this is the entire third word of the line, as a double
            normalization = Double.parseDouble(words[2]);
            break;
          case "B053F08":
            // this is the normalization frequency
            // once again the entire third word of the line as double
            normalFreq = Double.parseDouble(words[2]);
            break;
          case "B053F09":
            // the number of zeros listed in reponse pole/zero lines
            // again, this is the entire third word, as an int
            int numZero = Integer.parseInt(words[2]);
            zerosArr = new Complex[numZero];
            break;
          case "B053F14":
            // same as above line but for the number of poles
            int numPole = Integer.parseInt(words[2]);
            polesArr = new Complex[numPole];
            break;
          case "B053F10-13":
            // these are the lists of response zeros, in order with index,
            // real (double), imaginary (double), & corresponding error terms
            parseAndApplyComplex(line,zerosArr);
            break;
          case "B053F15-18":
            // as above but for poles
            parseAndApplyComplex(line, polesArr);
            break;
          case "B058F03":
            // gain stage sequence number; again, full third word as int
            // this is used to map the gain value to an index
            gainStage = Integer.parseInt(words[2]);
            break;
          case "B058F04":
            // gain value -- has to have had stage sequence already set
            // again, it's the third full word, this time as a double
            gain[gainStage] = Double.parseDouble(words[2]);
            // now reset gainStage to prevent it getting modified again
            gainStage = -1;
            break;
          }
          
          line = br.readLine();
        } // else
        
      } // end of file-read loop (EOF reached, line is null)
      
      zeros = Arrays.asList(zerosArr);
      poles = Arrays.asList(polesArr);
      
      br.close();
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
  }
  
  private void parseAndApplyComplex(String line, Complex[] array) {
    // reparse the line. why are we doing this?
    // if a number is negative, only one space between it and prev. number
    // and the previous split operation assumed > 2 spaces between numbers
    String[] words = line.split("[ ]{1,}");


    // index 0 is the identifier for the field types
    // index 1 is where in the list this zero or pole is
    // index 2 is the real part, and index 3 the imaginary
    // indices 4 and 5 are error terms (ignored)    
    int index = Integer.parseInt(words[1]);
    double realPart = Double.parseDouble(words[2]);
    double imagPart = Double.parseDouble(words[3]);
    array[index] = new Complex(realPart, imagPart);
  }
  
  
  
}


