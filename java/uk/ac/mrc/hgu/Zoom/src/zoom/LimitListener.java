package zoom;

import zoom.*;
import java.util.EventObject;

/**
 * Event interface for Zoom.
 * @author Nick Burton
 * @see LimitEvent
 * @see Zoom
 */
public interface LimitListener {

/**
 * @param LimitEvent event
 * @return void
 */
   public void  limitChanged(LimitEvent event) ;

} // interface LimitListener
