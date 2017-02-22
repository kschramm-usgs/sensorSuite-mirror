package asl.sensor;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.MultivariateJacobianFunction;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.RealVectorFormat;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.Pair;
import org.jfree.data.xy.XYSeries;

public class OrthogonalExperiment extends Experiment {

  double[] diffs;
  
  public OrthogonalExperiment() {
    super();

  }

  @Override
  protected void backend(final DataStore ds, final boolean freqSpace) {
    // TODO auto generated method stub
    
    long interval = ds.getXthLoadedBlock(1).getInterval();
    
    // assume the first two are the reference and the second two are the test?
    
    // TODO: replace with calls to get Xth block in case we expand plots again
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
        start(new double[] {0.,0.,0.,0.}).
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
    
    RealVector solution = optimum.getPoint();
    double sum = solution.getEntry(0) + solution.getEntry(1);
    solution = solution.mapDivide( sum );
    
    RealVectorFormat rvf = new RealVectorFormat("[","]",", ");
    System.out.println(rvf.format(solution));
    
    diffs = solution.toArray();
    
    double timeAtPoint = 0.;
    double tick = interval / TimeSeriesUtils.ONE_HZ_INTERVAL;
    double scale1 = diffs[0];
    double scale2 = diffs[1];
    
    XYSeries diffSrs = new XYSeries("Diff(" + testName + ", " + refName + ")");
    XYSeries diffRotSrs = new XYSeries("Diff(" + testName + ", Rotated Ref.)");
    
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
  

}
