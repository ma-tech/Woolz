package wsetter;

import wsetter.*;
import java.util.EventObject;

/**
 * Event interface for WSetter.
 * @author Nick Burton
 * @see LimitEvent
 * @see WSetter
 */
public interface LimitListener {

   public void  limitChanged(LimitEvent event) ;

} // interface LimitListener
