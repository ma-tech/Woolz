package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;


/**
 *   An entry in the AnatKey.
 *   KeyEntry has controls to allow the colour to be changed
 *   and the component visibility to be toggled off / on.
 */
public class KeyEntry extends KeyEntryGUI{

   protected Color _col = null;

   protected String _txt = "";

   protected int _indx = 0;

   protected boolean _visible = true;
   protected boolean _visible3D = true;
   protected boolean _removed = false;

   /**
    *   Constructs a 2D KeyEntry.
    */
   public KeyEntry(int indx) {
      this(indx, false);
   }
   
   /**
    *   Constructs a 2D or 3D KeyEntry depending
    *   upon the state of 'is3D'.
    */
   public KeyEntry(int indx, boolean is3D) {

      super(is3D);

      ButtonHandler0 buttonH0 = new ButtonHandler0();
      colBtn0.addActionListener(buttonH0);

      ButtonHandler1 buttonH1 = new ButtonHandler1();
      btn00.addActionListener(buttonH1);
      btn01.addActionListener(buttonH1);

      ButtonHandler2 buttonH2 = new ButtonHandler2();
      btn02.addActionListener(buttonH2);

      if(is3D) {
	 ButtonHandler3 buttonH3 = new ButtonHandler3();
	 btn03.addActionListener(buttonH3);
      }

      ButtonHandler4 buttonH4 = new ButtonHandler4();
      btn04.addActionListener(buttonH4);

      _indx = indx;

   }
   
//-------------------------------------------------------------
   /**
    *   Sets the anatomy component name for the KeyEntry.
    *   @param str the full name of the anatomy component.
    */
   public void setText(String str) {

      _txt = str;
      tf0.setText(_txt);
   }

//-------------------------------------------------------------
   /**
    *   Returns the name of the anatomy component represented by
    *   this EntryKey.
    *   @return _txt.
    */
   protected String getTxt(){
     return _txt;
   }

//-------------------------------------------------------------
   /**
    *   Removes the anatomy component from the KeyEntry.
    *   Changes the 'zap' button icon to a bent arrow.
    */
   public void setEntryRemoved() {

      tf0.setText("");
      tf0.setBackground(_notVizCol);
      btn02.setIcon(replaceIcon);
   }

//-------------------------------------------------------------
   /**
    *   Sets the visibility of the anatomy component.
    *   @param str the full name of the anatomy component.
    */
   public void setEntryVisible(boolean viz) {

      tf0.setText(getTxt());
      if(viz) {
	 tf0.setBackground(_vizCol);
      } else {
         tf0.setBackground(_notVizCol);
      }
   }

//-------------------------------------------------------------
   /**
    *   Sets the colour of the anatomy component
    *   to that returned by a standard color chooser dialog.
    *   The color chooser dialog is opened by clicking on the 
    *   Colour chooser button in the KeyEntry.
    *   @param col the new colour of the anatomy component.
    */
   public void setCol(Color col) {

      _col = col;
      colBtn0.setBackground(_col);
   }

//-------------------------------------------------------------
   /**
    *   Returns the colour of an anatomy component at
    *   the given index in the Anatomy Key.
    *   @param index the index.
    *   @return a new Color object.
    */
   public Color getCol(){
     return _col;
   }

//-------------------------------------------------------------
   /**
    *   Resets the KeyEntry to its empty state.
    */
   public void reset() {

      tf0.setText("");
      tf0.setBackground(_notVizCol);
      btn02.setIcon(null);
   }

//-------------------------------------------------------------
   /**
    *   Updates the KeyEntry with new / modified data.
    *   Array entries may be null.
    *   @param arr an array of AnatomyElements.
    *   @see sectionViewer.AnatomyElement
    */
/*
   public void update(AnatomyElement el) {

      if(el != null) {
	 if(el.isRemoved()) {
	    setEntryRemoved();
	 } else {
	    if(el.isVisible()) {
	       setEntryVisible(el.getDescription(), true);
	    } else {
	       setEntryVisible(el.getDescription(), false);
	    }
	 }
      }
   }
*/

//-------------------------------------------------------------
   /**
    *   Returns the width of a KeyEntry.
    *   @return _entryW
    */
   public static int getW() {
      return _entryW;
   }

//-------------------------------------------------------------
   /**
    *   Returns the height of a KeyEntry.
    *   @return _entryH
    */
   public static int getH() {
      return _entryH;
   }

//-------------------------------------------------------------
   /**
    *   Returns the unique index of a KeyEntry.
    *   @return _indx
    */
   public int getIndx() {
      return _indx;
   }

//-------------------------------------------------------------
   /**
    *   True if the entry has been removed.
    *   @return _removed
    */
   public boolean isRemoved() {
      return _removed;
   }

//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
   /**
    *   Event handler for the colour chooser button.
    *   This concatenates "colour" with the index number
    *   of the component as the parameter for a 'fireAction' command.
    *   An application  will implement an ActionListener to listen for
    *   these events and manage the updates required on the KeyEntry.
    */
   public class ButtonHandler0 implements ActionListener {

      public ButtonHandler0() {
      }

      public void actionPerformed(ActionEvent e) {
         if(e.getSource() == colBtn0) {
	    fireAction(new String("colour" + _indx));
         }
      }

   } // inner class ButtonHandler0

//-------------------------------------------------------------
   /**
    *   Event handler for the arrow (scroll) buttons for the textfield
    *   These scroll the anatomy name in the text field
    *   fully left or right.
    */
   public class ButtonHandler1 implements ActionListener {

      /*
         not sure how to calculate width of text in pixels
         as I can't seem to get the font from the JTextField
         hence END = a big number !!
      */
      final int START = 0;
      final int END = 4000;

      public ButtonHandler1() {
      }

      public void actionPerformed(ActionEvent e) {
         if(e.getSource() == btn00) {
	    tf0.setScrollOffset(START); // got to start of text
         } else if(e.getSource() == btn01) {
            tf0.setScrollOffset(END); // go to end of text
         }
      }

   } // inner class ButtonHandler1

//-------------------------------------------------------------
   /**
    *   Event handler for the 2D visibility button.
    *   This concatenates "makeVisible" (or makeInvisible) with
    *   the index number of the component as the
    *   parameter for a 'fireAction' command.
    *   An application  will implement an ActionListener to listen for
    *   these events and manage the updates required on the KeyEntry.
    */
   public class ButtonHandler2 implements ActionListener {

      public ButtonHandler2() {
      }

      public void actionPerformed(ActionEvent e) {
	 if(e.getSource() == btn02) {
	    if(!_visible) {
	       if(vizIcon != null) {
	            btn02.setIcon(vizIcon);
               }
	       fireAction(new String("makeVisible" + _indx));
	    } else {
	       if(notVizIcon != null) {
	            btn02.setIcon(notVizIcon);
               }
	       fireAction(new String("makeInvisible" + _indx));
	    }
	    _visible = !_visible;
         }
      }

   } // inner class ButtonHandler2

//-------------------------------------------------------------
   /**
    *   Event handler for the 3D visibility button.
    *   This concatenates "makeVisible" (or makeInvisible) with
    *   the index number of the component as the
    *   parameter for a 'fireAction' command.
    *   An application  will implement an ActionListener to listen for
    *   these events and manage the updates required on the KeyEntry.
    */
   public class ButtonHandler3 implements ActionListener {

      public ButtonHandler3() {
      }

      public void actionPerformed(ActionEvent e) {
	 if(e.getSource() == btn03) {
	    if(!_visible3D) {
	       if(viz3DIcon != null) {
	            btn03.setIcon(viz3DIcon);
               }
	       fireAction(new String("showThreeDAnat" + _indx));
	    } else {
	       if(notViz3DIcon != null) {
	            btn03.setIcon(notViz3DIcon);
               }
	       fireAction(new String("hideThreeDAnat" + _indx));
	    }
	    _visible3D = !_visible3D;
         }
      }

   } // inner class ButtonHandler3

//-------------------------------------------------------------
   /**
    *   Event handler for the zap button.
    *   This concatenates "zap" with the relevant index (row) number
    *   of the component as the parameter for a 'fireAction' command.
    *   An application  will implement an ActionListener to listen for
    *   these events and manage the updates required on the KeyEntry.
    */
   public class ButtonHandler4 implements ActionListener {

      public ButtonHandler4() {
      }

      public void actionPerformed(ActionEvent e) {
         if(e.getSource() == btn04) {
	    fireAction(new String("zap" + _indx));
         }
      }

   } // inner class ButtonHandler3

//-------------------------------------------------------------
// handle all objects that are interested in changes
//-------------------------------------------------------------

   // keep track of all the listeners to this model
   /**
    *   A list of ActionListeners which are
    *   listening for events fired from the KeyEntry.
    */
   protected EventListenerList actionListeners =
                             new EventListenerList();

  // add a listener to the register
   /**
    *   Registers an ActionListener
    *   with the EventListenerList.
    *   @param x an Event handler implementing ActionListener
    */
  public void addActionListener(ActionListener x) {
    actionListeners.add (ActionListener.class, x);
  }


  // remove a listener from the register
   /**
    *   Removes a previously registered ActionListener
    *   from the EventListenerList
    *   @param x an Event handler implementing ActionListener
    */
  public void removeActionListener(ActionListener x) {
    actionListeners.remove (ActionListener.class, x);
  }

   /**
    *   Fires an ActionEvent with the given Action Command.
    *   This allows an event handler to choose the appropriate action.
    */
   protected void fireAction(String cmd) {
   // Create the event:
   ActionEvent ae = new ActionEvent(this,
                                    ActionEvent.ACTION_PERFORMED,
				    cmd);
   // Get the listener list
   Object[] listeners =
     actionListeners.getListenerList();
   // Process the listeners last to first
   // List is in pairs, Class and instance
   for (int i
     = listeners.length-2; i >= 0; i -= 2) {
     if (listeners[i] == ActionListener.class) {
        ActionListener al = (ActionListener)listeners[i+1];
        al.actionPerformed(ae);
     }
   }
  } // fireAction

//-------------------------------------------------------------

} // class KeyEntry
