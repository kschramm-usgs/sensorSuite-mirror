package asl.sensor;

/**
 * Used to specify the type of transfer function transform for response
 * calculation, based on the associated parameters in the response file
 * @author akearns
 *
 */
public enum TransferFunction {

  /**
   * Specifies that response is calculated with 2*pi*f*i for frequency f
   * (Equivalent to Fourier transform given positive time values)
   */
  LAPLACIAN (2*Math.PI),
  /**
   * Specifies that response is calculated with f*i for frequency f
   */
  LINEAR (1.0);
  
  TransferFunction(double function) {
    this.function = function;
  }
  
  private double function;
  
  public double getFunction(){
    return function;
  }
  
}
