package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.net.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;

/**
 *   The GUI for KeyEntry.
 */
public class KeyEntryGUI extends JPanel{


   /** The system file separator ("/" or "\"). */
   private String SLASH = System.getProperty("file.separator");

   /** The colour chooser button for row 0 of the KeyEntry   */
   protected JButton colBtn0 = null;

   /** The left arrow button for row 0 of the KeyEntry   */
   protected JButton btn00 = null;

   /** The right arrow button for row 0 of the KeyEntry   */
   protected JButton btn01 = null;

   /** The scrollable text field for row 0 of the KeyEntry   */
   protected ScrollableTextField tf0 = null;

   /** The visibility button of a KeyEntry   */
   protected JButton btn02 = null;

   /** The zap button of a KeyEntry   */
   protected JButton btn03 = null;

   /** The 3D visibility button of a KeyEntry   */
   protected JButton btn04 = null;

   /** The 'left arrow' icon for the textfield scroll button  */
   protected ImageIcon leftIcon = null;
   /** The 'right arrow' icon for the textfield scroll button  */
   protected ImageIcon rightIcon = null;
   /** The 'viz' icon (the 2D anatomy component is visible)  */
   protected ImageIcon vizIcon = null;
   /** The 'notViz' icon (the 2D anatomy component is not visible)  */
   protected ImageIcon notVizIcon = null;
   /** The 'viz3D' icon (the 3D anatomy component is visible)  */
   protected ImageIcon viz3DIcon = null;
   /** The 'notViz3D' icon (the 3D anatomy component is not visible)  */
   protected ImageIcon notViz3DIcon = null;
   /** The 'zap' icon (toggles between zap & replace)  */
   protected ImageIcon zapIcon = null;
   /** The 'replace' icon (toggles between zap & replace)  */
   protected ImageIcon replaceIcon = null;


   /** The tool tip for the colour chooser button */
   protected String colTipText = "choose colour";
   /** The tool tip for the left arrow button */
   protected String leftTipText = "scroll to start of text";
   /** The tool tip for the right arrow button */
   protected String rightTipText = "scroll to end of text";
   /** The tool tip for the visibility button */
   protected String vizTip2Text = "toggle 2D visibility";
   /** The tool tip for the visibility button */
   protected String vizTip3Text = "toggle 3D visibility";
   /** The tool tip for the zap button */
   protected String zapTipText = "remove from key";

   /**
    *   The background colour of the textfield for the anatomy component name
    *   when it is currently visible. (Visible means that it will be
    *   displayed if the section intersects the component).
    */

   protected Color _vizCol;
   /**
    *   The background colour of the textfield for the anatomy component name
    *   when it is currently not visible. (Visible means that it will not be
    *   displayed even if the section intersects the component).
    */
   protected Color _notVizCol;

   protected boolean _is3D;
   protected static final int _entryW = 200;
   //protected static final int _entryH = 35;
   protected static final int _entryH = 25;
//-------------------------------------------------------------
   /**
    *   Constructs a KeyEntryGUI.
    */
   protected KeyEntryGUI() {
      _is3D = false;
      makeGUI();
   }

   /**
    *   Constructs a KeyEntryGUI with 3D visibility control if
    *   is3D is true.
    *   @param is3D true if 3D visibility control is to be added.
    */
   protected KeyEntryGUI(boolean is3D) {
      _is3D = is3D;
      makeGUI();
   }

//-------------------------------------------------------------
   /**
    *   Called by the constructor. All of the work of building the GUI
    *   is done here.
    */
   private void makeGUI() {

      int btnW1 = 18; // left, right button panels
      int btnW2 = 20; // colour, viz & zap button panels
      int btnW3 = 25; // combined left & right button panel
      int gap = 1;
      int W1 = 2*btnW1+btnW3+2*gap;

      try {
        String imgPath = "images/";
        URL urlIcon = KeyEntryGUI.class.getResource(imgPath + "left.gif");
        leftIcon = new ImageIcon(urlIcon);
        urlIcon = KeyEntryGUI.class.getResource(imgPath + "right.gif");
        rightIcon = new ImageIcon(urlIcon);
        urlIcon = KeyEntryGUI.class.getResource(imgPath + "viz.gif");
        vizIcon = new ImageIcon(urlIcon);
        urlIcon = KeyEntryGUI.class.getResource(imgPath + "notViz.gif");
        notVizIcon = new ImageIcon(urlIcon);
        urlIcon = KeyEntryGUI.class.getResource(imgPath + "viz.gif");
        viz3DIcon = new ImageIcon(urlIcon);
        urlIcon = KeyEntryGUI.class.getResource(imgPath + "notViz.gif");
        notViz3DIcon = new ImageIcon(urlIcon);
        urlIcon = KeyEntryGUI.class.getResource(imgPath + "zap.gif");
        zapIcon = new ImageIcon(urlIcon);
        urlIcon = KeyEntryGUI.class.getResource(imgPath + "replace.gif");
        replaceIcon = new ImageIcon(urlIcon);
      }
      catch (Exception ex) {
      }

      _vizCol = new Color(250,250,250);
      _notVizCol = new Color(200,200,200);

//......................................
      this.setLayout(new BorderLayout());
      this.setPreferredSize(new Dimension(0, _entryH));
      //this.setBackground(Color.green);
//......................................
      /**
       *   Outer container for 1 row of the anatomy key
       *   which maintains the preferred size.
       *   <br>Contains kPanel0.
       */
      JPanel kPanel = new JPanel();

      /**
       *   Inner container for 1 row of the anatomy key.
       *   <br>Contains kPanel00, txtPanel0, kPanel02.
       */
      JPanel kPanel0 = new JPanel();

      /**
       *   Container for <em>colour chooser</em> and <em>scroll</em> buttons.
       *   <br>Contains colPanel0, navPanel0.
       *   <br>Added to kPanel0.
       */
      JPanel kPanel00 = new JPanel();

      /**
       *   Container for <em>visibility</em> controls
       *   <br>Contains
       *   <br>Added to kPanel02.
       */
      JPanel kPanel01 = new JPanel();

      /**
       *   Container for <em>visibility</em> and <em>zap</em> controls
       *   <br>Contains kPanel02, btnPanel03.
       *   <br>Added to kPanel0.
       */
      JPanel kPanel02 = new JPanel();

      kPanel.setLayout(new BorderLayout(gap,gap));
      //kPanel.setBackground(Color.orange);
      kPanel.setPreferredSize(new Dimension(_entryW,_entryH));

      kPanel0.setLayout(new BorderLayout(gap,gap));
      //kPanel0.setBackground(Color.pink);
      kPanel0.setPreferredSize(new Dimension(0,_entryH));

      kPanel00.setLayout(new BorderLayout(gap,gap));
      //kPanel00.setBackground(Color.blue);
      kPanel00.setPreferredSize(new Dimension(W1,0));
      kPanel00.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      kPanel01.setLayout(new BorderLayout(gap,gap));
      /*
      kPanel01.setBackground(Color.magenta);
      kPanel01.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */

      kPanel02.setLayout(new BorderLayout(gap,gap));
      //kPanel02.setBackground(Color.magenta);
      kPanel02.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      colBtn0 = new JButton();
      JPanel colPanel0 = new JPanel();
      colPanel0.setLayout(new BorderLayout(gap,gap));
      //colPanel0.setPreferredSize(new Dimension(btnW2,0));
      colPanel0.setPreferredSize(new Dimension(btnW3,0));
      colPanel0.add(colBtn0);
      JPanel txtPanel0 = new JPanel();
      txtPanel0.setLayout(new BorderLayout(gap,gap));
      txtPanel0.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      tf0 = new ScrollableTextField("", 50);
      tf0.setEditable(false);
      tf0.setBackground(_vizCol);
      txtPanel0.add(tf0, BorderLayout.CENTER);
      JPanel navPanel0 = new JPanel();
      navPanel0.setPreferredSize(new Dimension(2*btnW1,0));
      //navPanel0.setLayout(new BoxLayout(navPanel0, BoxLayout.X_AXIS));
      navPanel0.setLayout(new BorderLayout(gap,gap));
      JPanel btnPanel00 = new JPanel();
      btnPanel00.setPreferredSize(new Dimension(btnW1,0));
      btnPanel00.setLayout(new BorderLayout(gap,gap));
      btn00 = new JButton();
      if(leftIcon != null) {
	 btn00.setIcon(leftIcon);
      }
      btnPanel00.add(btn00, BorderLayout.CENTER);
      JPanel btnPanel01 = new JPanel();
      btnPanel01.setPreferredSize(new Dimension(btnW1,0));
      btnPanel01.setLayout(new BorderLayout(gap,gap));
      btn01 = new JButton();
      if(leftIcon != null) {
	 btn01.setIcon(rightIcon);
      }
      btnPanel01.add(btn01, BorderLayout.CENTER);
      navPanel0.add(btnPanel00, BorderLayout.WEST);
      navPanel0.add(btnPanel01, BorderLayout.EAST);

      JPanel btnPanel02 = new JPanel();
      btnPanel02.setPreferredSize(new Dimension(btnW2,0));
      btnPanel02.setLayout(new BorderLayout(gap,gap));
      btn02 = new JButton();
      if(vizIcon != null) {
	 btn02.setIcon(vizIcon);
      }
      btnPanel02.add(btn02, BorderLayout.CENTER);

      JPanel btnPanel03 = new JPanel();
      btnPanel03.setPreferredSize(new Dimension(btnW2,0));
      btnPanel03.setLayout(new BorderLayout(gap,gap));
      btn03 = new JButton();
      if(viz3DIcon != null) {
	 btn03.setIcon(viz3DIcon);
      }
      btnPanel03.add(btn03, BorderLayout.CENTER);

      JPanel btnPanel04 = new JPanel();
      btnPanel04.setPreferredSize(new Dimension(btnW2,0));
      btnPanel04.setLayout(new BorderLayout(gap,gap));
      btn04 = new JButton();
      if(zapIcon != null) {
	 btn04.setIcon(zapIcon);
      }
      btnPanel04.add(btn04, BorderLayout.CENTER);

      kPanel00.add(colPanel0, BorderLayout.WEST);
      kPanel00.add(navPanel0, BorderLayout.EAST);
      if(_is3D) {
	 kPanel01.add(btnPanel02, BorderLayout.WEST);
	 kPanel01.add(btnPanel03, BorderLayout.EAST);
      } else {
         kPanel01.setPreferredSize(new Dimension(2*btnW2+5*gap,0));
	 kPanel01.add(btnPanel02, BorderLayout.CENTER);
      }
      kPanel02.add(kPanel01, BorderLayout.WEST);
      kPanel02.add(btnPanel04, BorderLayout.EAST);
      kPanel0.add(kPanel00, BorderLayout.WEST);
      kPanel0.add(txtPanel0, BorderLayout.CENTER);
      kPanel0.add(kPanel02, BorderLayout.EAST);

      //kPanel.add(kPanel0, BorderLayout.NORTH); // set height to preferred size
      kPanel.add(kPanel0, BorderLayout.CENTER); // set height to expand
      this.add(kPanel, BorderLayout.CENTER); // allows text field to expand
//......................................
      setToolTips();
//......................................

   } // makeKeyEntryGUI()

//......................................

   /**
    *   Attaches tool tip text to the controls in the KeyEntry.
    */
   private void setToolTips() {

      ToolTipManager.sharedInstance().setInitialDelay( 1500 );
      ToolTipManager.sharedInstance().setDismissDelay( 1500 );

      colBtn0.setToolTipText(colTipText);
      btn00.setToolTipText(leftTipText);
      btn01.setToolTipText(rightTipText);
      btn02.setToolTipText(vizTip2Text);
      btn03.setToolTipText(vizTip3Text);
      btn04.setToolTipText(zapTipText);

   }

//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
   /**
    *   Abstract event handler for colour chooser button.
    */
   protected abstract class ButtonHandler0 implements ActionListener {
   }

   /**
    *   Abstract event handler for arrow (scroll) buttons.
    */
   protected abstract class ButtonHandler1 implements ActionListener {
   }

   /**
    *   Abstract event handler for visibility button.
    */
   protected abstract class ButtonHandler2 implements ActionListener {
   }

   /**
    *   Abstract event handler for visibility button.
    */
   protected abstract class ButtonHandler3 implements ActionListener {
   }

   /**
    *   Abstract event handler for zap button.
    */
   protected abstract class ButtonHandler4 implements ActionListener {
   }

} // class KeyEntryGUI
