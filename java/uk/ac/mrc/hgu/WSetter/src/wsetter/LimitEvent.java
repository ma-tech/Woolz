package wsetter;

import wsetter.*;
import java.util.EventObject;

/**
 * Event class for WSetter.
 * Created when WSetter model's limits are changed.
 * @author Nick Burton
 * @see WSetterModel WSetter
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
