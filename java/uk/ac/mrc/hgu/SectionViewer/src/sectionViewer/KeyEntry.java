package sectionViewer;

import java.awt.*;
import java.awt.event.*;
import javax.swing.event.*;


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
      colBtn.addActionListener(buttonH0);

      ButtonHandler1 buttonH1 = new ButtonHandler1();
      leftBtn.addActionListener(buttonH1);
      rightBtn.addActionListener(buttonH1);

      ButtonHandler2 buttonH2 = new ButtonHandler2();
      vis2DBtn.addActionListener(buttonH2);

      if(is3D) {
	 ButtonHandler3 buttonH3 = new ButtonHandler3();
	 vis3DBtn.addActionListener(buttonH3);
      }

      ButtonHandler4 buttonH4 = new ButtonHandler4();
      zapBtn.addActionListener(buttonH4);

      _indx = indx;

   }

//-------------------------------------------------------------
   /**
    *   Sets the anatomy component name for the KeyEntry.
    *   @param str the full name of the anatomy component.
    */
   public void setText(String str) {

      _txt = str;
      textF.setText(_txt);
   }

//-------------------------------------------------------------
   /**
    *   Returns the name of the anatomy component represented by
    *   this EntryKey.
    *   @return _txt.
    */
   public String getTxt(){
     return _txt;
   }

//-------------------------------------------------------------
   /**
    *   Sets the text for the anatomy component.
    */
   public void setEntryText() {
      textF.setText(getTxt());
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
      colBtn.setBackground(_col);
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

      textF.setText("");
      vis2DBtn.setIcon(null);
   }

//-------------------------------------------------------------
   /**
    *   Returns the width of a KeyEntry.
    *   @return _entryW
    */
   public static int getW() {
      return _entryW;
   }

//......................................
   /**
    *   Sets the width of a KeyEntry.
    */
   protected void setKeyWidth(int w) {
      this.setPreferredSize(new Dimension(w, _entryH));
      this.setMinimumSize(new Dimension(w, _entryH));
      this.setMaximumSize(new Dimension(10000, _entryH));
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
    *   Sets the icon for 3D visibility control.
    *   The KeyEntry will not fire an event.
    *   @param state true if icon is 'visible'.
    */
   protected void set3DVisIcon(boolean visible) {
      if(visible) {
         vis3DBtn.setIcon(viz3DIcon);
      } else {
         vis3DBtn.setIcon(notViz3DIcon);
      }
      _visible3D = visible;
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
         if(e.getSource() == colBtn) {
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
         if(e.getSource() == leftBtn) {
	    textF.setScrollOffset(START); // got to start of text
         } else if(e.getSource() == rightBtn) {
            textF.setScrollOffset(END); // go to end of text
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
	 if(e.getSource() == vis2DBtn) {
	    vis2DBtn.setToolTipText(null);
	    if(!_visible) {
	       if(hideIcon != null) {
	            vis2DBtn.setIcon(hideIcon);
		    vis2DBtn.setToolTipText(hide2TipText);
               }
	       fireAction(new String("makeVisible" + _indx));
	    } else {
	       if(showIcon != null) {
	            vis2DBtn.setIcon(showIcon);
		    vis2DBtn.setToolTipText(show2TipText);
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
	 if(e.getSource() == vis3DBtn) {
	    vis3DBtn.setToolTipText(null);
	    if(!_visible3D) {
	       if(hideIcon != null) {
	            vis3DBtn.setIcon(hideIcon);
		    vis3DBtn.setToolTipText(hide3TipText);
               }
	       fireAction(new String("showThreeDAnat" + _indx));
	    } else {
	       if(showIcon != null) {
	            vis3DBtn.setIcon(showIcon);
		    vis3DBtn.setToolTipText(show3TipText);
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
         if(e.getSource() == zapBtn) {
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
