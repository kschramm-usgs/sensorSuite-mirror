package asl.sensor;

import java.awt.Dimension;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.Scrollable;

public class VScrollPanel extends JPanel implements Scrollable {

  List<JPanel> subpanelList;
  
  public VScrollPanel() {
    super();
    subpanelList = new ArrayList<JPanel>();
  }
  
  public void add(JPanel jp) {
    super.add(jp);
    subpanelList.add(jp);
  }
  
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
    return false;
  }

  
  
}
