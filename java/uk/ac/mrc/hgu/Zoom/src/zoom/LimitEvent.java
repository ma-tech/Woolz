package zoom;

import zoom.*;
import java.util.EventObject;

/**
 * Event class for Zoom.
 * Created when Zoom model's limits are changed.
 * @author Nick Burton
 * @see ZoomModel Zoom
 */
public class LimitEvent extends EventObject {

//----------------------------------------------------
  /**
   * 
   */
   public LimitEvent(Object source) {
      super(source);
   }

} // class LimitEvent
