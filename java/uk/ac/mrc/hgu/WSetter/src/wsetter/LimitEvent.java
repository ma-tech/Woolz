package wsetter;

import wsetter.*;
import java.util.EventObject;

/**
 * Event class for WSetter bean.
 * @author Nick Burton
 */
public class LimitEvent extends EventObject {

//----------------------------------------------------
  /**
   * Constructor
   */
   public LimitEvent(Object source) {
      super(source);
   }

} // class LimitEvent
