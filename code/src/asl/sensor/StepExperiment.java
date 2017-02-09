package asl.sensor;

import java.util.List;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.leastsquares.EvaluationRmsChecker;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.
        leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.
        leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;

public class StepExperiment extends Experiment{

  double f, h; //corner and damping of output (uncorrected)
  double fCorr, hCorr; // fit parameters to turn output into cal input
  
  int trimLength;
  double[] freqs;
  Complex[] sensorFFTSeries; // FFT of step cal from sensor
  double[] stepCalSeries; // time series of raw step cal function
  
  final double STEP_FACTOR = Math.nextUp( Math.nextUp( (double) 1.0 ) );
  
  public StepExperiment() {
    super();
  }
  
  @Override
  protected void backend(final DataStore ds, final boolean freqSpace) {
    
    // assume that the first block is the raw step calibration
    // the raw calibration is defined as not having an associated response
    DataBlock stepCalRaw = ds.getXthLoadedBlock(1);
    
    stepCalSeries = new double[stepCalRaw.size()];
    

    
    // resulting inverse fft should be this length
    trimLength = stepCalRaw.size();
    long interval = stepCalRaw.getInterval();
    
    double[] stepCalUnfiltered = new double[stepCalRaw.size()];
    for (int i = 0; i < trimLength; ++i) {
      stepCalUnfiltered[i] = stepCalRaw.getData().get(i).doubleValue();
    }
    
    
    /*
    Complex[] fftStepCal = FFTResult.simpleFFT(stepCalRaw);
        double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / interval;
    
    double freqPerCell = sps / (fftStepCal.length/2);
    
    Complex[] firstHalfStepCal = new Complex[fftStepCal.length/2];
    
    for (int i = 0; i < firstHalfStepCal.length; ++i) {
      if ( i > 0.01 / sps * firstHalfStepCal.length) {
        firstHalfStepCal[i] = Complex.ZERO;
      } else {
        firstHalfStepCal[i] = fftStepCal[i];
      }
    }
    
    double[] stepCalFiltered = 
                            FFTResult.inverseFFT(firstHalfStepCal, trimLength);
    */
    
    // stepCalSeries = stepCalFiltered;
    
    stepCalSeries = TimeSeriesUtils.normalize(stepCalUnfiltered);
    
    // FFTResult.detrend(toDetrend);
    
    // stepCalSeries = TimeSeriesUtils.normalize(filteredStepCal);
    // stepCalSeries = filteredStepCal;
    
    // but we want the response and the data of the cal result
    int outIdx = ds.getXthFullyLoadedIndex(1);
    
    // if someone did load a raw cal with the response, then we wouldn't
    // get a different block with the second call above, so we get the 
    // next loaded block/response pair
    if ( ds.getBlock(outIdx).getName().equals( stepCalRaw.getName() ) ) {
      outIdx = ds.getXthFullyLoadedIndex(2);
    }

    // get data of the result of the step calibration
    DataBlock sensorOutput = ds.getBlock(outIdx);
    // long interval = sensorOutput.getInterval();
    InstrumentResponse ir = ds.getResponse(outIdx);
    Complex pole = ir.getPoles().get(0);
    
    f = 1. / (2 * Math.PI / pole.abs() ); // corner frequency
    h = Math.abs( pole.getReal() / pole.abs() ); // damping
    
    double[] params = new double[]{f, h};
    
    // get FFT of datablock timeseries, deconvolve with response
    FFTResult sensorsFFT = FFTResult.simpleSingleSidedFFT(sensorOutput);
    sensorFFTSeries = sensorsFFT.getFFT();
    freqs = sensorsFFT.getFreqs();
    
    double[] toPlot = calculate(params);

    toPlot = TimeSeriesUtils.normalize(toPlot);
    
    long start = 0L; // was db.startTime();
    long now = start;
    XYSeries xys = new XYSeries("STEP *^(-1) RESP");
    XYSeries scs = new XYSeries( stepCalRaw.getName() ); 
    for (double point : toPlot) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      xys.add(seconds, point);
      now += interval;
    }
    now = start;
    for (Number point : stepCalSeries ) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      scs.add(seconds, point);
      now += interval;
    }
    
    // next we'll want to find the parameters to fit the plots
    // to the inputted data
    xySeriesData.addSeries(scs);
    xySeriesData.addSeries(xys);
    
    // next step: curve fitting
    RealVector startVector = MatrixUtils.createRealVector(params);
    RealVector observedComponents = MatrixUtils.createRealVector(stepCalSeries);
    
    
    ConvergenceChecker<LeastSquaresProblem.Evaluation> svc = 
        new EvaluationRmsChecker(1E-8, 1E-8);
        
    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(startVector).
        target(observedComponents).
        model( getJacobianFunction() ).
        lazyEvaluation(false).
        maxEvaluations(1000000000).
        maxIterations(1000000000).
        checker(svc).
        build();
    
    LeastSquaresProblem.Evaluation initEval = lsp.evaluate(startVector);
    System.out.println("INITIAL GUESS RMS: " +  initEval.getRMS() );
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer();
        //withCostRelativeTolerance(1.0E-17).
        //withParameterRelativeTolerance(1.0E-16);
    
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    
    System.out.println("FIT PARAMS RMS: " +  optimum.getRMS() );
    
    double[] newParams = optimum.getPoint().toArray();
    
    fCorr = newParams[0];
    hCorr = newParams[1];
    
    double[] fitPlot = FFTResult.demean( calculate(newParams) );
    fitPlot = TimeSeriesUtils.normalize(fitPlot);
    start = 0L; // was db.startTime();
    now = start;
    XYSeries bfs = new XYSeries("BEST FIT PLOT");
    for (double point : fitPlot) {
      double seconds = (double) now / TimeSeriesUtils.ONE_HZ_INTERVAL;
      bfs.add(seconds, point);
      now += interval;
    }
    
    xySeriesData.addSeries(bfs);
    
    
  }
  
  public double[] calculate(double[] beta) {
    
    double f = beta[0];
    double h = beta[1];
    
    // use these terms to calculate the poles
    
    // term inside the square root in the calculations of p1, p2
    // (h^2-1)
    Complex tempResult = new Complex( Math.pow(h, 2) - 1 );
    
    double omega = 2 * Math.PI * f; // omega_0
    
    // casting h as a complex number to enable these calculations
    Complex hCast = new Complex(h);
    
    // - (h + sqrt(h^2-1))
    Complex pole1 = hCast.add( tempResult.sqrt() ).multiply(-1);
    pole1 = pole1.multiply(omega); 
      // cannot modify complex numbers "in-place"; they MUST be (re)assigned
    
    // - (h - sqrt(h^2-1))
    Complex pole2 = hCast.subtract( tempResult.sqrt() ).multiply(-1);
    pole2 = pole2.multiply(omega);
    
    // calculate the FFT of the response
    Complex[] respFFT = new Complex[sensorFFTSeries.length]; // array of resps
    double max = 0.0;
    // don't let denominator be zero
    respFFT[0] = Complex.ONE;
    for (int i = 1; i < respFFT.length; ++i) {

      // replaced freqs[i] with 
      // 2*pi*i*f
      Complex factor = new Complex(0, 2*Math.PI*freqs[i]); 
      
      // (2*pi*i*f - p1) * (2*pi*f*i - p2)
      Complex denom = factor.subtract(pole1).multiply( factor.subtract(pole2) );
      
      respFFT[i] = factor.divide(denom);
      
      
      if (respFFT[i].abs() > max) {
        max = respFFT[i].abs();
      }
      
    }
    
    Complex[] toDeconvolve = new Complex[sensorFFTSeries.length];
    
    // deconvolving response is dividing fft(signal) by fft(response)
    
    for (int i = 0; i < sensorFFTSeries.length; ++i) {
      // the conjugate of the response, used twice in deconvolution
      Complex conjResp = respFFT[i].conjugate();
      double aboveZero = 0.008*max;  // term to keep the denominator above 0
      
      // resp * conj(resp) + 0.008(max(|resp|))
      Complex denom = respFFT[i].multiply(conjResp).add(aboveZero);
      
      // deconvolving the response from output
      // fft * conj(resp) / (resp * conjResp)+0.008(max(|resp|))
      toDeconvolve[i] = sensorFFTSeries[i].multiply(conjResp).divide(denom);
    }
    
    double[] returnValue =  FFTResult.inverseFFT(toDeconvolve, trimLength);
    
    return returnValue;
    
  }
  
  private Pair<RealVector, RealMatrix> jacobian(RealVector variables) {
    
    // approximate through forward differences
    double[][] jacobian = new double[trimLength][2];
    
    double f1 = variables.getEntry(0);
    double h1 = variables.getEntry(1);
    double f2 = f1 * STEP_FACTOR;
    double h2 = h1 * STEP_FACTOR;
    
    double[] fInit = FFTResult.demean( calculate(new double[]{f1, h1}) );
    double[] diffOnF = FFTResult.demean( calculate(new double[]{f2, h1}) );
    double[] diffOnH = FFTResult.demean( calculate(new double[]{f1, h2}) );
    
    fInit = TimeSeriesUtils.normalize(fInit);
    diffOnF = TimeSeriesUtils.normalize(diffOnF);
    diffOnH = TimeSeriesUtils.normalize(diffOnH);
    
    for (int i = 0; i < trimLength; ++i) {
      jacobian[i][0] = (diffOnF[i] - fInit[i]) / (f2 - f1);
      jacobian[i][1] = (diffOnH[i] - fInit[i]) / (h2 - h1);
    }
    
    RealMatrix jMat = MatrixUtils.createRealMatrix(jacobian);
    RealVector fnc = MatrixUtils.createRealVector(fInit);
    
    return new Pair<RealVector, RealMatrix>(fnc, jMat);
  }
  
  public MultivariateJacobianFunction getJacobianFunction() {
    return new MultivariateJacobianFunction() {
      private static final long serialVersionUID = -8673650298627399464L;
      public Pair<RealVector, RealMatrix> value(RealVector point) {
          return jacobian(point);
      }
    };
  }
  
  public double[] value(double[] point) {
    return calculate(point);
  }
  
  public double[] getFitCornerAndDamping() {
    return new double[]{fCorr, hCorr};
  }

  public double[] getCornerAndDamping() {
    return new double[]{f, h};
  }
  
}
