package asl.sensor;

import java.awt.Dimension;
import java.awt.Rectangle;

import javax.swing.JPanel;
import javax.swing.Scrollable;

public class VScrollPanel extends JPanel implements Scrollable {

  
  /**
   * 
   */
  private static final long serialVersionUID = -4358405572297138013L;

  boolean scaleVerticallyNotScroll = false;
  
  @Override
  public Dimension getPreferredScrollableViewportSize() {
    return getPreferredSize();
  }

  @Override
  public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation,
      int direction) {
    return 1;
  }

  @Override
  public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation,
      int direction) {
    return 1;
  }

  @Override
  public boolean getScrollableTracksViewportWidth() {
    return true;
  }

  @Override
  public boolean getScrollableTracksViewportHeight() {
    return scaleVerticallyNotScroll;
  }

  public void setScrollableTracksViewportHeight(boolean bool) {
    scaleVerticallyNotScroll = bool;
  }
  
}