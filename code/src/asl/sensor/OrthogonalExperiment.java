package asl.sensor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import 
org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import 
org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.RealVectorFormat;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;

public class OrthogonalExperiment extends Experiment {

  double[] diffs;
  double angle;
  
  public OrthogonalExperiment() {
    super();

  }

  @Override
  protected void backend(final DataStore ds, final boolean freqSpace) {
    // TODO auto generated method stub
    
    long interval = ds.getXthLoadedBlock(1).getInterval();
    
    // assume the first two are the reference and the second two are the test?
    
    // we just need four timeseries, don't actually care about response
    DataBlock refLH1Block = ds.getXthLoadedBlock(1);
    DataBlock refLH2Block = ds.getXthLoadedBlock(2);
    DataBlock testLH1Block = ds.getXthLoadedBlock(3);
    DataBlock testLH2Block = ds.getXthLoadedBlock(4);    
    
    List<Number> refLH1 = new ArrayList<Number>( refLH1Block.getData() );
    String refName = refLH1Block.getName();
    List<Number> refLH2 = new ArrayList<Number>( refLH2Block.getData() );
    List<Number> testLH1 = new ArrayList<Number>( testLH1Block.getData() );
    String testName = testLH1Block.getName();
    List<Number> testLH2 = new ArrayList<Number>( testLH2Block.getData() );
    
    FFTResult.detrend(refLH1);
    FFTResult.detrend(refLH2);
    FFTResult.detrend(testLH1);
    FFTResult.detrend(testLH2);
    
    refLH1 = TimeSeriesUtils.normalize(refLH1);
    refLH2 = TimeSeriesUtils.normalize(refLH2);
    testLH1 = TimeSeriesUtils.normalize(testLH1);
    testLH2 = TimeSeriesUtils.normalize(testLH2);
    
    double sps = TimeSeriesUtils.ONE_HZ_INTERVAL / refLH1Block.getInterval();
    double low = 1./8;
    double high = 1./4;
    
    refLH1 = FFTResult.bandFilter(refLH1, sps, low, high);
    refLH2 = FFTResult.bandFilter(refLH2, sps, low, high);
    testLH1 = FFTResult.bandFilter(testLH1, sps, low, high);
    testLH2 = FFTResult.bandFilter(testLH2, sps, low, high);
    
    // System.out.println(refLH1);
    
    // will assume all vectors are trimmed to the same size
    int len = refLH1.size();
    double[] rhsArray = new double[len * 2];
    double[][] lhsArray = new double[len * 2][4];
    
    
    for (int i = 0; i < len; ++i) {
      rhsArray[i] = refLH1.get(i).doubleValue();
      rhsArray[i + len] = refLH2.get(i).doubleValue();
      
      lhsArray[i][0] = testLH1.get(i).doubleValue();
      lhsArray[i][1] = testLH2.get(i).doubleValue();
      lhsArray[i + len][2] = testLH1.get(i).doubleValue();
      lhsArray[i + len][3] = testLH2.get(i).doubleValue();
    }
    
    RealVector rhsVector = MatrixUtils.createRealVector(rhsArray);
    RealMatrix lhsMatrix = MatrixUtils.createRealMatrix(lhsArray);
    
    MultivariateJacobianFunction mjf = new MultivariateJacobianFunction() {
      public Pair<RealVector, RealMatrix> value(final RealVector point) {
        
        RealVector result = lhsMatrix.operate(point);
        return new Pair<RealVector, RealMatrix>(result, lhsMatrix);
      }
    };
    
    LeastSquaresProblem lsp = new LeastSquaresBuilder().
        start(new double[] {1.,0.,0.,1.}).
        model(mjf).
        target(rhsVector).
        maxEvaluations(Integer.MAX_VALUE).
        maxIterations(Integer.MAX_VALUE).
        lazyEvaluation(false).
        build();
    
    LeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer().
        withCostRelativeTolerance(1.0E-15).
        withParameterRelativeTolerance(1.0E-15);
    
    LeastSquaresOptimizer.Optimum optimum = optimizer.optimize(lsp);
    
    RealVector fitSolution = optimum.getPoint();
    // double divisor = solution.getEntry(0) + solution.getEntry(1);
        // solution.getMaxValue();
    // solution = solution.mapDivide(divisor);
    
    RealVectorFormat rvf = new RealVectorFormat("[","]",", ");
    System.out.println(rvf.format(fitSolution));
    
    diffs = fitSolution.toArray();
    
    double timeAtPoint = 0.;
    double tick = interval / TimeSeriesUtils.ONE_HZ_INTERVAL;
    double scale1 = diffs[0];
    double scale2 = diffs[1];
    
    XYSeries diffSrs = new XYSeries("Diff(" + testName + ", " + refName + ")");
    XYSeries diffRotSrs = new XYSeries("Diff(" + testName + ", Rotated Ref.)");
    
    double[] angs = new double[4];
    
    angs[0] = FastMath.toDegrees( FastMath.acos(diffs[0]) ) % 360;
    angs[1] = FastMath.toDegrees( FastMath.asin(diffs[1]) ) % 360;
    angs[2] = FastMath.toDegrees( -FastMath.asin(diffs[2]) ) % 360;
    angs[3] = FastMath.toDegrees( FastMath.acos(diffs[3]) ) % 360;
    
    
    double avg = 0.;
    for (int i = 0; i < angs.length; ++i) {
      if ( angs[i] < 0 ) {
        angs[i] += 360; // it's not modulo in java; results can be negative
      }
      
      avg += angs[i];
    }
    
    avg /= angs.length;
    
    angle = avg;
    
    double sDev = 0.;
    
    for (double ang : angs) {
      sDev += Math.pow( avg - ang, 2 );
    }
    
    sDev = Math.sqrt(sDev / 4);
    
    System.out.println(Arrays.toString(angs));
    System.out.println(sDev);
    
    for (int i = 0; i < len; ++i) {
      double rotSubt = rhsArray[i] * scale1 + rhsArray[len + i] * scale2;
      
      double unrotatedDiff = lhsArray[i][0] - rhsArray[i];
      double rotatedDiff = lhsArray[i][0] - rotSubt;
      
      diffSrs.add(timeAtPoint, unrotatedDiff);
      diffRotSrs.add(timeAtPoint, rotatedDiff);
      
      timeAtPoint += tick;
    }
    
    xySeriesData.addSeries(diffSrs);
    xySeriesData.addSeries(diffRotSrs);
    
  }

  public double[] getSolutionParams() {
    return diffs;
  }
  
  public double getFitAngle() {
    return angle;
  }
  

}
