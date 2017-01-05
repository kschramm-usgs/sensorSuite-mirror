package asl.sensor.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import asl.sensor.FFTResult;

public class FFTResultTest {

  @Test
  public void rangeCopyTest() {
  
    Number[] numbers = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    
    int low = 5;
    int high = 9;
    
    List<Number> numList = Arrays.asList(numbers);
    List<Number> subseq = new ArrayList<Number>(numList.subList(low, high));
    
    for (int i = 0; i < subseq.size(); ++i) {
      int fullListIdx = i + low;
      assertEquals( numList.get(fullListIdx), subseq.get(i) );
    }
    
    for (int i = 0; i < subseq.size(); ++i) {
      Number temp = subseq.get(i);
      temp = 2000;
      assertNotEquals(subseq.get(i),temp); // can't try to change "in-place"
      subseq.set(i, 100);
    }
    
    for (int i = 0; i < subseq.size(); ++i) {
      int fullListIdx = i + low;
      assertNotEquals( numList.get(fullListIdx), subseq.get(i) );
    }
    
  }
  
  @Test
  public void demeaningTest() {
    
    // tests that demean does what it says it does and that
    // the results are applied in-place
    
    Number[] numbers = {1,2,3,4,5};
    
    List<Number> numList = Arrays.asList(numbers);
    List<Number> demeaned = new ArrayList<Number>(numList);
    
    FFTResult.demean(demeaned);
    
    for (int i = 0; i < numList.size(); ++i) {
      assertEquals(demeaned.get(i), numList.get(i).doubleValue()-3);
    }
    
  }
  
  @Test
  public void detrendingLinearTest() {
    
    Number[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
    
    List<Number> toDetrend = Arrays.asList(x);
    FFTResult.detrend(toDetrend);
    
    for (Number num : toDetrend) {
      assertEquals(num.doubleValue(), 0.0, 0.001);
    }
    
  }

  
  @Test
  public void detrendingCycleTest() {
    
    Number[] x = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
        18, 19, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 
        3, 2, 1 };
    
    List<Number> toDetrend = Arrays.asList(x);
    
    Number[] answer = { -9d, -8d, -7d, -6d, -5d, -4d, -3d, -2d, -1d, 0d, 1d, 2d,
        3d, 4d, 5d, 6d, 7d, 8d, 9d, 10d, 9d, 8d, 7d, 6d, 5d, 4d, 3d, 2d, 1d, 0d,
        -1d, -2d, -3d, -4d, -5d, -6d, -7d, -8d, -9d };

    
    FFTResult.detrend(toDetrend);
    
    for (int i = 0; i < x.length; i++) {
      assertEquals(
          new Double(Math.round(x[i].doubleValue())), 
          new Double(answer[i].doubleValue()));
    }
    
  }
  
  @Test
  public final void cosineTaperTest() throws Exception {
    Number[] x = { 5, 5, 5, 5, 5 };
    List<Number> toTaper = Arrays.asList(x);
    Double[] tapered = { 0d, 4.5d, 5d, 4.5d, 0d };

    double power = FFTResult.cosineTaper(toTaper, 0.25);
    
    assertEquals(new Double(Math.round(power)), new Double(4));
    
    for (int i = 0; i < x.length; i++) {
      // round to the first decimal (multiply by 10, round, divide by 10)
      Double result = 
          new Double(Math.round(toTaper.get(i).doubleValue()*10d)/10d);
      assertEquals(result, tapered[i]);
    }
  }
  
}
