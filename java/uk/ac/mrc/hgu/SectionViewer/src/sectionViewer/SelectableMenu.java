package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

public class SelectableMenu extends JMenu {

   ActionEvent event;
   String modifiersStr;

   public SelectableMenu() {
      super();
   }
   public SelectableMenu(String str) {
      super(str);
   }

   public void processMouseEvent (MouseEvent e,
                                  MenuElement[] path,
				  MenuSelectionManager man) {

      // comparing bit patterns doesn't seem to work
      // Shift is 64 but modifiers returns 17
      modifiersStr = e.getMouseModifiersText(e.getModifiers());

      switch (e.getButton()) {

      case MouseEvent.BUTTON1:
	 if(modifiersStr.indexOf("Shift") != -1) {
	    event = new ActionEvent(this,
				    ActionEvent.ACTION_PERFORMED,
				    this.getActionCommand());
	    man.clearSelectedPath();
	    fireEvent();
	 }
	 break;
      case MouseEvent.BUTTON2:
      case MouseEvent.BUTTON3:
	 event = new ActionEvent(this,
				 ActionEvent.ACTION_PERFORMED,
				 this.getActionCommand());
	 man.clearSelectedPath();
	 fireEvent();
	 break;
      default:
	 break;
      }
   }

//-------------------------------------------------------------
//-------------------------------------------------------------
   // keep track of all the listeners to this 'model'
   protected EventListenerList actionListeners =
                             new EventListenerList();

//-------------------------------------------------------------
  // add a listener to the register
  public void addActionListener(ActionListener x) {
    actionListeners.add (ActionListener.class, x);
  }

//-------------------------------------------------------------
  // remove a listener from the register
  public void removeActionListener(ActionListener x) {
    actionListeners.remove (ActionListener.class, x);
  }

//-------------------------------------------------------------
   protected void fireEvent() {
   // Get the listener list
   Object[] listeners =
     actionListeners.getListenerList();
   // Process the listeners last to first
   // List is in pairs, Class and instance
   for (int i
     = listeners.length-2; i >= 0; i -= 2) {
     if (listeners[i] == ActionListener.class) {
        ActionListener al = (ActionListener)listeners[i+1];
        al.actionPerformed(event);
     }
   }
  } // fireEvent

}

 
