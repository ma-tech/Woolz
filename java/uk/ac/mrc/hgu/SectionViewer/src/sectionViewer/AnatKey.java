package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import java.io.*;

import uk.ac.mrc.hgu.Wlz.*;

public class AnatKey extends AnatKeyGUI{

   static protected  AnatKey _instance = null;
   static protected  final int _nrows = 6;
//-------------------------------------------------------------
   // private constructor for singleton pattern
   private AnatKey(String str) {
      super(str);

      ButtonHandler1 buttonH1 = new ButtonHandler1();
      btn00.addActionListener(buttonH1);
      btn10.addActionListener(buttonH1);
      btn20.addActionListener(buttonH1);
      btn30.addActionListener(buttonH1);
      btn40.addActionListener(buttonH1);
      btn50.addActionListener(buttonH1);
      btn01.addActionListener(buttonH1);
      btn11.addActionListener(buttonH1);
      btn21.addActionListener(buttonH1);
      btn31.addActionListener(buttonH1);
      btn41.addActionListener(buttonH1);
      btn51.addActionListener(buttonH1);

      ButtonHandler2 buttonH2 = new ButtonHandler2();
      btn02.addActionListener(buttonH2);
      btn12.addActionListener(buttonH2);
      btn22.addActionListener(buttonH2);
      btn32.addActionListener(buttonH2);
      btn42.addActionListener(buttonH2);
      btn52.addActionListener(buttonH2);

      ButtonHandler3 buttonH3 = new ButtonHandler3();
      colBtn0.addActionListener(buttonH3);
      colBtn1.addActionListener(buttonH3);
      colBtn2.addActionListener(buttonH3);
      colBtn3.addActionListener(buttonH3);
      colBtn4.addActionListener(buttonH3);
      colBtn5.addActionListener(buttonH3);

      ChBoxHandler cbxH = new ChBoxHandler();
      cbx0.addItemListener(cbxH);
      cbx1.addItemListener(cbxH);
      cbx2.addItemListener(cbxH);
      cbx3.addItemListener(cbxH);
      cbx4.addItemListener(cbxH);
      cbx5.addItemListener(cbxH);

   }
   
//-------------------------------------------------------------
   // makes sure only 1 instance of the class is created
   static public AnatKey instance() {
      if(_instance == null) {
         _instance = new AnatKey("anatomy key");
      }
      return _instance;
   }
//-------------------------------------------------------------
   public int getNRows() {
      return _nrows;
   }
//-------------------------------------------------------------
   public static void setText(String str, int indx) {

      switch (indx) {
	 case 0:
	    tf0.setText(str);
	    break;
	 case 1:
	    tf1.setText(str);
	    break;
	 case 2:
	    tf2.setText(str);
	    break;
	 case 3:
	    tf3.setText(str);
	    break;
	 case 4:
	    tf4.setText(str);
	    break;
	 case 5:
	    tf5.setText(str);
	    break;
	 default:
	    break;
      } // switch
   }

//-------------------------------------------------------------
   public static void setText(String[] strArr) {
      tf0.setText(strArr[0]);
      tf1.setText(strArr[1]);
      tf2.setText(strArr[2]);
      tf3.setText(strArr[3]);
      tf4.setText(strArr[4]);
      tf5.setText(strArr[5]);
   }

//-------------------------------------------------------------
   public static void setEntryRemoved(int indx, String str) {

      switch(indx) {
      case 0:
	 tf0.setText("");
	 tf0.setBackground(notvis);
	 btn02.setIcon(replaceIcon);
         break;
      case 1:
	 tf1.setText("");
	 tf1.setBackground(notvis);
	 btn12.setIcon(replaceIcon);
         break;
      case 2:
	 tf2.setText("");
	 tf2.setBackground(notvis);
	 btn22.setIcon(replaceIcon);
         break;
      case 3:
	 tf3.setText("");
	 tf3.setBackground(notvis);
	 btn32.setIcon(replaceIcon);
         break;
      case 4:
	 tf4.setText("");
	 tf4.setBackground(notvis);
	 btn42.setIcon(replaceIcon);
         break;
      case 5:
	 tf5.setText("");
	 tf5.setBackground(notvis);
	 btn52.setIcon(replaceIcon);
         break;
      default:
         break;
      }
   }

//-------------------------------------------------------------
   public static void setEntryVisible(int indx, String str) {

      switch(indx) {
      case 0:
         tf0.setText(str);
         tf0.setBackground(vis);
	 cbx0.setSelected(true);
	 btn02.setIcon(zapIcon);
         break;
      case 1:
         tf1.setText(str);
         tf1.setBackground(vis);
	 cbx1.setSelected(true);
	 btn12.setIcon(zapIcon);
         break;
      case 2:
         tf2.setText(str);
         tf2.setBackground(vis);
	 cbx2.setSelected(true);
	 btn22.setIcon(zapIcon);
         break;
      case 3:
         tf3.setText(str);
         tf3.setBackground(vis);
	 cbx3.setSelected(true);
	 btn32.setIcon(zapIcon);
         break;
      case 4:
         tf4.setText(str);
         tf4.setBackground(vis);
	 cbx4.setSelected(true);
	 btn42.setIcon(zapIcon);
         break;
      case 5:
         tf5.setText(str);
         tf5.setBackground(vis);
	 cbx5.setSelected(true);
	 btn52.setIcon(zapIcon);
         break;
      default:
         break;
      }
   }

//-------------------------------------------------------------
   public static void setEntryInvisible(int indx, String str) {

      switch(indx) {
      case 0:
         tf0.setText(str);
         tf0.setBackground(notvis);
	 cbx0.setSelected(false);
	 btn02.setIcon(zapIcon);
         break;
      case 1:
         tf1.setText(str);
         tf1.setBackground(notvis);
	 cbx1.setSelected(false);
	 btn12.setIcon(zapIcon);
         break;
      case 2:
         tf2.setText(str);
         tf2.setBackground(notvis);
	 cbx2.setSelected(false);
	 btn22.setIcon(zapIcon);
         break;
      case 3:
         tf3.setText(str);
         tf3.setBackground(notvis);
	 cbx3.setSelected(false);
	 btn32.setIcon(zapIcon);
         break;
      case 4:
         tf4.setText(str);
         tf4.setBackground(notvis);
	 cbx4.setSelected(false);
	 btn42.setIcon(zapIcon);
         break;
      case 5:
         tf5.setText(str);
         tf5.setBackground(notvis);
	 cbx5.setSelected(false);
	 btn52.setIcon(zapIcon);
         break;
      default:
         break;
      }
   }

//-------------------------------------------------------------
   public static void setCol(Color col, int indx) {

      _cols[indx] = col.getBlue() |
                    col.getGreen() << 8 |
		    col.getRed() << 16 |
		    col.getAlpha() << 24;

      switch(indx) {
      case 0:
         colBtn0.setBackground(new Color(_cols[indx]));
	 break;
      case 1:
         colBtn1.setBackground(new Color(_cols[indx]));
	 break;
      case 2:
         colBtn2.setBackground(new Color(_cols[indx]));
	 break;
      case 3:
         colBtn3.setBackground(new Color(_cols[indx]));
	 break;
      case 4:
         colBtn4.setBackground(new Color(_cols[indx]));
	 break;
      case 5:
         colBtn5.setBackground(new Color(_cols[indx]));
	 break;
      default:
         break;
      } // switch
   }

//-------------------------------------------------------------
   public static void reset() {

      tf0.setText("");
      tf0.setBackground(notvis);
      cbx0.setSelected(false);
      btn02.setIcon(null);

      tf1.setText("");
      tf1.setBackground(notvis);
      cbx1.setSelected(false);
      btn12.setIcon(null);

      tf2.setText("");
      tf2.setBackground(notvis);
      cbx2.setSelected(false);
      btn22.setIcon(null);

      tf3.setText("");
      tf3.setBackground(notvis);
      cbx3.setSelected(false);
      btn32.setIcon(null);

      tf4.setText("");
      tf4.setBackground(notvis);
      cbx4.setSelected(false);
      btn42.setIcon(null);

      tf5.setText("");
      tf5.setBackground(notvis);
      cbx5.setSelected(false);
      btn52.setIcon(null);
   }

//-------------------------------------------------------------
   public static void update(AnatomyElement[] arr) {

      AnatomyElement el = null;

      for(int i=0; i<_nrows; i++) {
         el = arr[i];
         if(el != null) {
	    if(el.isRemoved()) {
	       setEntryRemoved(i, el.getDescription());
	    } else {
	       if(el.isVisible()) {
		  setEntryVisible(i, el.getDescription());	    
	       } else {
		  setEntryInvisible(i, el.getDescription());	    
	       }
	    }
	 }
      } // for
   }

//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
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
         } else if(e.getSource() == btn10) {
            tf1.setScrollOffset(START); 
         } else if(e.getSource() == btn11) {
            tf1.setScrollOffset(END);
         } else if(e.getSource() == btn20) {
            tf2.setScrollOffset(START); 
         } else if(e.getSource() == btn21) {
            tf2.setScrollOffset(END);
         } else if(e.getSource() == btn30) {
            tf3.setScrollOffset(START); 
         } else if(e.getSource() == btn31) {
            tf3.setScrollOffset(END);
         } else if(e.getSource() == btn40) {
            tf4.setScrollOffset(START); 
         } else if(e.getSource() == btn41) {
            tf4.setScrollOffset(END);
         } else if(e.getSource() == btn50) {
            tf5.setScrollOffset(START); 
         } else if(e.getSource() == btn51) {
            tf5.setScrollOffset(END);
         }
      }

   } // inner class ButtonHandler1

//-------------------------------------------------------------
   public class ButtonHandler2 implements ActionListener {

      String cmd1 = "zap";

      public ButtonHandler2() {
      }

      public void actionPerformed(ActionEvent e) {
         if(e.getSource() == btn02) {
	    fireAction(new String(cmd1+"0"));
         } else if(e.getSource() == btn12) {
	    fireAction(new String(cmd1+"1"));
         } else if(e.getSource() == btn22) {
	    fireAction(new String(cmd1+"2"));
         } else if(e.getSource() == btn32) {
	    fireAction(new String(cmd1+"3"));
         } else if(e.getSource() == btn42) {
	    fireAction(new String(cmd1+"4"));
         } else if(e.getSource() == btn52) {
	    fireAction(new String(cmd1+"5"));
         }
      }

   } // inner class ButtonHandler2

//-------------------------------------------------------------
   public class ButtonHandler3 implements ActionListener {

      String cmd1 = "colour";

      public ButtonHandler3() {
      }

      public void actionPerformed(ActionEvent e) {
         if(e.getSource() == colBtn0) {
	    fireAction(new String(cmd1+"0"));
         } else if(e.getSource() == colBtn1) {
	    fireAction(new String(cmd1+"1"));
         } else if(e.getSource() == colBtn2) {
	    fireAction(new String(cmd1+"2"));
         } else if(e.getSource() == colBtn3) {
	    fireAction(new String(cmd1+"3"));
         } else if(e.getSource() == colBtn4) {
	    fireAction(new String(cmd1+"4"));
         } else if(e.getSource() == colBtn5) {
	    fireAction(new String(cmd1+"5"));
         }
      }

   } // inner class ButtonHandler3

//-------------------------------------------------------------
   public class ChBoxHandler implements ItemListener {

      String cmd1 = "makeVisible";
      String cmd2 = "makeInvisible";

      public ChBoxHandler() {
      }

      public void itemStateChanged(ItemEvent e) {
	 if(e.getSource() == cbx0) {
	    if(e.getStateChange() == ItemEvent.SELECTED) {
	       fireAction(new String(cmd1+"0"));
	    } else {
	       fireAction(new String(cmd2+"0"));
	    }
         } else if(e.getSource() == cbx1) {
	    if(e.getStateChange() == ItemEvent.SELECTED) {
	       fireAction(new String(cmd1+"1"));
	    } else {
	       fireAction(new String(cmd2+"1"));
	    }
         } else if(e.getSource() == cbx2) {
	    if(e.getStateChange() == ItemEvent.SELECTED) {
	       fireAction(new String(cmd1+"2"));
	    } else {
	       fireAction(new String(cmd2+"2"));
	    }
         } else if(e.getSource() == cbx3) {
	    if(e.getStateChange() == ItemEvent.SELECTED) {
	       fireAction(new String(cmd1+"3"));
	    } else {
	       fireAction(new String(cmd2+"3"));
	    }
         } else if(e.getSource() == cbx4) {
	    if(e.getStateChange() == ItemEvent.SELECTED) {
	       fireAction(new String(cmd1+"4"));
	    } else {
	       fireAction(new String(cmd2+"4"));
	    }
         } else if(e.getSource() == cbx5) {
	    if(e.getStateChange() == ItemEvent.SELECTED) {
	       fireAction(new String(cmd1+"5"));
	    } else {
	       fireAction(new String(cmd2+"5"));
	    }
         }
      }

   } // inner class chBoxHandler

//-------------------------------------------------------------
// handle all objects that are interested in changes
//-------------------------------------------------------------

   // keep track of all the listeners to this model
   protected EventListenerList actionListeners =
                             new EventListenerList();

  // add a listener to the register
  public void addActionListener(ActionListener x) {
    actionListeners.add (ActionListener.class, x);
  }


  // remove a listener from the register
  public void removeActionListener(ActionListener x) {
    actionListeners.remove (ActionListener.class, x);
  }

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

  public static void main(String argv[]) {

	AnatKey _key;

        _key = AnatKey.instance();
	_key.setSize(300,500);
	_key.pack();
	_key.setVisible(true);

  }

} // class AnatKey
