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
  LAPLACIAN,
  /**
   * Specifies that response is calculated with f*i for frequency f
   */
  LINEAR;
  
}
