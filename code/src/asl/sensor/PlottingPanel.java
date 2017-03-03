package asl.sensor;

import java.awt.Dimension;
import java.awt.Rectangle;

import javax.swing.Scrollable;

public class PlottingPanel implements Scrollable {

  @Override
  public Dimension getPreferredScrollableViewportSize() {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation,
      int direction) {
    // TODO Auto-generated method stub
    return 0;
  }

  @Override
  public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation,
      int direction) {
    // TODO Auto-generated method stub
    return 0;
  }

  @Override
  public boolean getScrollableTracksViewportWidth() {
    // TODO Auto-generated method stub
    return false;
  }

  @Override
  public boolean getScrollableTracksViewportHeight() {
    // TODO Auto-generated method stub
    return false;
  }

}
