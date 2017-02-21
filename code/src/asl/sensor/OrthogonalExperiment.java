package asl.sensor;

import java.util.List;

import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.RealVectorFormat;
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
    List<Number> refLH1 = ds.getBlock(0).getData();
    String refName = ds.getBlock(0).getName();
    List<Number> refLH2 = ds.getBlock(1).getData();
    List<Number> testLH1 = ds.getBlock(2).getData();
    String testName = ds.getBlock(2).getName();
    List<Number> testLH2 = ds.getBlock(3).getData();
    
    // will assume all vectors are trimmed to the same size
    int len = refLH1.size();
    double[] rhsArray = new double[len * 2];
    double[][] lhsArray = new double[len * 2][4];
    
    for (int i = 0; i < len * 2; ++i) {
      rhsArray[i] = refLH1.get(i).doubleValue();
      rhsArray[i + len] = refLH2.get(i).doubleValue();
      
      lhsArray[i][0] = testLH1.get(i).doubleValue();
      lhsArray[i][1] = testLH2.get(i).doubleValue();
      lhsArray[i + len][2] = testLH1.get(i).doubleValue();
      lhsArray[i + len][3] = testLH2.get(i).doubleValue();
    }
    
    RealVector rhsVector = MatrixUtils.createRealVector(rhsArray);
    RealMatrix lhsMatrix = MatrixUtils.createRealMatrix(lhsArray);
    
    // QR decomposition solver uses least-squares to fit variables
    DecompositionSolver solver = new QRDecomposition(lhsMatrix).getSolver();
    RealVector solution = solver.solve(rhsVector);
    // scale the solution
    double maxVal = solution.getMaxValue();
    solution = solution.mapDivide(maxVal);
    RealVectorFormat rvf = new RealVectorFormat("[","]",",");
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
