package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;

/**
 *   Dialogue with fields for x, y, z coordinates.
 */
public class PointEntry extends PointEntryGUI{
//public class PointEntry extends JFrame{

  private double[] _values;
  private boolean damper = true;

//-------------------------------------------------------------
  /**
   *   Creates a PointEntry with the given Title.
   *   @param str the title of the PointEntry Dialogue.
   */
  public PointEntry(String str) {
    super(str);

    ButtonHandler BH = new ButtonHandler();
    _okButton.addActionListener(BH);
    _cancelButton.addActionListener(BH);
    TextFieldHandler TFH = new TextFieldHandler();
    tfX.addActionListener(TFH);
    tfY.addActionListener(TFH);
    tfZ.addActionListener(TFH);
    _values = new double[3];
    _values[0] = 0.0;
    _values[1] = 0.0;
    _values[2] = 0.0;
  }
//-------------------------------------------------------------
  /**
   *   Returns the x, y, z coordinates.
   *   @return the x, y, z coordinates.
   */
  public double[] getValues() {
    return _values;
  }
//-------------------------------------------------------------
  /**
   *   Returns the PointEntry.
   *   @return a reference to the PointEntry.
   */
  public PointEntry getThis() {
    return this;
  }

//-------------------------------------------------------------
  /**
   *   Initialises the PointEntry with the given x, y, z coordinates.
   *   @param vals an array containing x, y, z coordinates.
   */
  public void setInitVals(double vals[]) {
    /*
       find out where the decimal point is and set
       length of string for 1 place of decimals
    */
    int pntpos = 0;

    String xstr = Double.toString(vals[0]);
    String ystr = Double.toString(vals[1]);
    String zstr = Double.toString(vals[2]);

    pntpos = xstr.indexOf(".");
    if((pntpos != -1) && (xstr.length() >= pntpos+2)) {
       tfX.setText(xstr.substring(0,pntpos+2));
    } else {
       tfX.setText(xstr);
    }

    pntpos = ystr.indexOf(".");
    if((pntpos != -1) && (ystr.length() >= pntpos+2)) {
       tfY.setText(ystr.substring(0,pntpos+2));
    } else {
       tfY.setText(ystr);
    }

    pntpos = zstr.indexOf(".");
    if((pntpos != -1) && (zstr.length() >= pntpos+2)) {
       tfZ.setText(zstr.substring(0,pntpos+2));
    } else {
       tfZ.setText(zstr);
    }
  }

//-------------------------------------------------------------
  /**
   *   Event handler for PointEntry dialogue.
   *   Called if the user presses <return> inside
   *   one of the PointEntry text fields.
   */
  class TextFieldHandler implements ActionListener {

     public TextFieldHandler() {
     }

     public void actionPerformed(ActionEvent e) {

	if((e.getSource().equals(tfX)) ||
	   (e.getSource().equals(tfY)) ||
	   (e.getSource().equals(tfZ)))  {

	   _okButton.doClick();
	}

     }
  }
//-------------------------------------------------------------
  /**
   *   Event handler for PointEntry dialogue.
   *   Called if the user presses the 'OK' or 'Cancel'
   *   button in the PointEntry.
   */
  class ButtonHandler implements ActionListener {

    String xstr = "";
    String ystr = "";
    String zstr = "";

    public ButtonHandler() {
    }

    public void actionPerformed(ActionEvent e) {

      if(e.getSource().equals(_okButton)) {
        // do validity check here
        xstr = tfX.getText();
        ystr = tfY.getText();
        zstr = tfZ.getText();

        _values[0] = (Double.valueOf(xstr)).doubleValue();
        _values[1] = (Double.valueOf(ystr)).doubleValue();
        _values[2] = (Double.valueOf(zstr)).doubleValue();

        getThis().setVisible(false);
        getThis().dispose();
        fireChange();

      } else if(e.getSource().equals(_cancelButton)) {
        getThis().setVisible(false);
        getThis().dispose();
      }
    }
  }

//-------------------------------------------------------------
//-------------------------------------------------------------
// handle all _objects that are interested in changes
//-------------------------------------------------------------
  // keep track of all the listeners to this 'model'
  /**
   *   A list of ChangeListeners which are
   *   listening for events fired from the PointEntry.
   */
  protected EventListenerList changeListeners =
      new EventListenerList();

//-------------------------------------------------------------
  // add a listener to the register
  /**
   *   Adds a ChangeListener to the EventListenerList.
   */
  public void addChangeListener(ChangeListener x) {
    changeListeners.add (ChangeListener.class, x);

    // bring it up to date with current state
    x.stateChanged(new ChangeEvent(this));
  }

//-------------------------------------------------------------
  // remove a listener from the register
  /**
   *   Removes a ChangeListener from the EventListenerList.
   */
  public void removeChangeListener(ChangeListener x) {
    changeListeners.remove (ChangeListener.class, x);
  }

//-------------------------------------------------------------
  /**
   *   Fires a ChangeEvent from the PointEntry.
   */
  protected void fireChange() {
    // Create the event:
    ChangeEvent ce = new ChangeEvent(this);
    // Get the listener list
    Object[] listeners =
        changeListeners.getListenerList();
    // Process the listeners last to first
    // List is in pairs, Class and instance
    for (int i
         = listeners.length-2; i >= 0; i -= 2) {
      if (listeners[i] == ChangeListener.class) {
        ChangeListener cl = (ChangeListener)listeners[i+1];
        cl.stateChanged(ce);
      }
    }
  } // fireChange

//======================================================

} // class PointEntry
