package zoom;
import zoom.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

public class NBButton extends JButton {

  protected EventListenerList myMouseListeners = null;

   /* constructors */
  public NBButton() {
    super();
    myMouseListeners = new EventListenerList();
    if(myMouseListeners == null) {
      System.out.println("myMouseListeners null in constructor");
    }

  }
  public NBButton(Icon icon) {
    super(icon);
    myMouseListeners = new EventListenerList();
    if(myMouseListeners == null) {
      System.out.println("myMouseListeners null in constructor");
    }
  }
  public NBButton(String str) {
    super(str);
    myMouseListeners = new EventListenerList();
    if(myMouseListeners == null) {
      System.out.println("myMouseListeners null in constructor");
    }
  }

  public void processMouseEvent (MouseEvent e) {
     System.out.println("NBButton processing MouseEvent");
     fireEvent(e);
     fireMyEvent(e);
  }

//-------------------------------------------------------------
   // keep track of all my listeners to this 'model'
//-------------------------------------------------------------
  // add a listener to the register
  public void addMyMouseListener(MouseListener x) {
    if(myMouseListeners == null) {
       System.out.println("event list was null");
       return;
    } else {
       System.out.println("adding listener");
       myMouseListeners.add (MouseListener.class, x);
       System.out.println("now "+
              myMouseListeners.getListenerCount()+" listeners");
    }
  }

//-------------------------------------------------------------
  // remove a listener from the register
  public void removeMyMouseListener(MouseListener x) {
    if(myMouseListeners == null) {
       System.out.println("event list was null");
       return;
    } else {
       System.out.println("removing listener");
       myMouseListeners.remove (MouseListener.class, x);
       System.out.println("now "+
              myMouseListeners.getListenerCount()+" listeners");
    }
  }

//-------------------------------------------------------------
  protected void fireEvent(MouseEvent e) {
     // Process the listeners last to first
     // List is in pairs, Class and instance
     MouseListener[] mouseListeners = this.getMouseListeners();
     if(mouseListeners == null) {
        System.out.println("no MouseListeners");
	return;
     } else {
        System.out.println(mouseListeners.length+" MouseListeners");
     }
     for (int i
	   = mouseListeners.length-1; i >= 0; i--) {
	MouseListener ml = mouseListeners[i];
	switch (e.getID()) {
	   case MouseEvent.MOUSE_PRESSED:
	      ml.mousePressed(e);
	      break;
	   case MouseEvent.MOUSE_RELEASED:
	      ml.mouseReleased(e);
	      break;
	   case MouseEvent.MOUSE_ENTERED:
	      ml.mouseEntered(e);
	      break;
	   case MouseEvent.MOUSE_EXITED:
	      ml.mouseExited(e);
	      break;
	   case MouseEvent.MOUSE_CLICKED:
	      ml.mouseClicked(e);
	      break;
	   default:
	      break;
	} // switch
     }
  } // fireEvent


//-------------------------------------------------------------
   protected void fireMyEvent(MouseEvent e) {
   // Get the listener list
   Object[] listeners =
     myMouseListeners.getListenerList();
   // Process the listeners last to first
   // List is in pairs, Class and instance
     if(listeners == null) {
        System.out.println("no myMouseListeners");
	return;
     } else {
        System.out.println((listeners.length)/2+" myMouseListeners");
     }
   for (int i
     = listeners.length-2; i >= 0; i -= 2) {
     if (listeners[i] == MouseListener.class) {
        MouseListener ml = (MouseListener)listeners[i+1];
	switch (e.getID()) {
	   case MouseEvent.MOUSE_PRESSED:
	      ml.mousePressed(e);
	      break;
	   case MouseEvent.MOUSE_RELEASED:
	      ml.mouseReleased(e);
	      break;
	   case MouseEvent.MOUSE_ENTERED:
	      ml.mouseEntered(e);
	      break;
	   case MouseEvent.MOUSE_EXITED:
	      ml.mouseExited(e);
	      break;
	   case MouseEvent.MOUSE_CLICKED:
	      ml.mouseClicked(e);
	      break;
	   default:
	      break;
	} // switch
     }
   }
  } // fireMyEvent

  public int totalListeners() {
     return myMouseListeners.getListenerCount();
  }

}

