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
   protected JButton colBtn = null;

   /** The left arrow button for row 0 of the KeyEntry   */
   protected JButton leftBtn = null;

   /** The right arrow button for row 0 of the KeyEntry   */
   protected JButton rightBtn = null;

   /** The scrollable text field for row 0 of the KeyEntry   */
   protected ScrollableTextField textF = null;

   /** The visibility button of a KeyEntry   */
   protected JButton vis2DBtn = null;

   /** The zap button of a KeyEntry   */
   protected JButton zapBtn = null;

   /** The 3D visibility button of a KeyEntry   */
   protected JButton vis3DBtn = null;

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
   /** The 'zap' icon (deletes an entry from the key)  */
   protected ImageIcon zapIcon = null;


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
   protected String zapTipText = "delete from key";

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
      this(false);
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
      int w = 0;


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
       *   <br>Contains innerPanel.
       */
      JPanel outerPanel = new JPanel();
      outerPanel.setLayout(new BorderLayout(gap,gap));
      outerPanel.setPreferredSize(new Dimension(_entryW,_entryH));

      /**
       *   Inner container for 1 row of the anatomy key.
       *   <br>Contains nameScrollPanel, buttonPanel.
       */
      JPanel innerPanel = new JPanel();
      innerPanel.setLayout(new BorderLayout(gap,gap));
      innerPanel.setPreferredSize(new Dimension(0,_entryH));
      //innerPanel.setBackground(Color.pink);

      //................................
      /**
       *   Contains leftBtnPanel, namePanel, rightBtnPanel..
       */
      JPanel scrollNamePanel = new JPanel();
      scrollNamePanel.setLayout(new BorderLayout(gap,gap));
      scrollNamePanel.setPreferredSize(new Dimension(0,_entryH));
      scrollNamePanel.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      /**
       *   Contains leftBtn.
       */
      JPanel leftBtnPanel = new JPanel();
      leftBtnPanel.setPreferredSize(new Dimension(btnW1,0));
      leftBtnPanel.setLayout(new BorderLayout(gap,gap));

      /**
       *   Contains rightBtn.
       */
      JPanel rightBtnPanel = new JPanel();
      rightBtnPanel.setPreferredSize(new Dimension(btnW1,0));
      rightBtnPanel.setLayout(new BorderLayout(gap,gap));

      /**
       *   Contains textF.
       */
      JPanel namePanel = new JPanel();
      namePanel.setLayout(new BorderLayout(gap,gap));
      namePanel.setPreferredSize(new Dimension(0,_entryH));

      //................................
      /**
       *   Contains vis2DBtn.
       */
      JPanel vis2DPanel = new JPanel();
      vis2DPanel.setLayout(new BorderLayout(gap,gap));
      vis2DPanel.setPreferredSize(new Dimension(btnW3,0));

      /**
       *   Contains vis3DBtn.
       */
      JPanel vis3DPanel = new JPanel();
      vis3DPanel.setLayout(new BorderLayout(gap,gap));
      vis3DPanel.setPreferredSize(new Dimension(btnW3,0));
      
      /**
       *   Contains vis2DPanel, vis3DPanel.
       */
      JPanel vis2D3DPanel = new JPanel();
      vis2D3DPanel.setLayout(new BorderLayout(gap,gap));
      //vis2D3DPanel.setBackground(Color.yellow);
      w = vis2DPanel.getPreferredSize().width;
      if(_is3D) {
         w += vis3DPanel.getPreferredSize().width;
	 vis2D3DPanel.setPreferredSize(new Dimension(w,0));
      } else {
	 vis2D3DPanel.setPreferredSize(new Dimension(w,0));
      }

      //................................
      /**
       *   Contains colBtn.
       */
      JPanel colPanel = new JPanel();
      colPanel.setLayout(new BorderLayout(gap,gap));
      colPanel.setPreferredSize(new Dimension(btnW3,0));

      /**
       *   Contains zapBtn.
       */
      JPanel zapPanel = new JPanel();
      zapPanel.setLayout(new BorderLayout(gap,gap));
      zapPanel.setPreferredSize(new Dimension(btnW3,0));

      /**
       *   Contains colPanel, zapPanel.
       */
      JPanel colZapPanel = new JPanel();
      colZapPanel.setLayout(new BorderLayout(gap,gap));
      w = colPanel.getPreferredSize().width;
      w += zapPanel.getPreferredSize().width;
      colZapPanel.setPreferredSize(new Dimension(w,0));


      //................................
      /**
       *   Contains vis2D3DPanel, colZapPanel.
       */
      JPanel buttonPanel = new JPanel();
      buttonPanel.setLayout(new BorderLayout(gap,gap));
      //buttonPanel.setBackground(Color.cyan);
      w = vis2D3DPanel.getPreferredSize().width;
      w += colZapPanel.getPreferredSize().width;
      buttonPanel.setPreferredSize(new Dimension(w,_entryH));
      //................................
//-----------------------------------------------------

      colBtn = new JButton();

      textF = new ScrollableTextField("", 50);
      textF.setEditable(false);
      textF.setBackground(_vizCol);
      
      leftBtn = new JButton();
      if(leftIcon != null) {
	 leftBtn.setIcon(leftIcon);
      }
      
      rightBtn = new JButton();
      if(leftIcon != null) {
	 rightBtn.setIcon(rightIcon);
      }

      vis2DBtn = new JButton();
      if(vizIcon != null) {
	 vis2DBtn.setIcon(vizIcon);
      }

      vis3DBtn = new JButton();
      if(viz3DIcon != null) {
	 vis3DBtn.setIcon(viz3DIcon);
      }

      zapBtn = new JButton();
      if(zapIcon != null) {
	 zapBtn.setIcon(zapIcon);
      }
//-----------------------------------------------------
      leftBtnPanel.add(leftBtn, BorderLayout.CENTER);
      namePanel.add(textF, BorderLayout.CENTER);
      rightBtnPanel.add(rightBtn, BorderLayout.CENTER);

      scrollNamePanel.add(leftBtnPanel, BorderLayout.WEST);
      scrollNamePanel.add(namePanel, BorderLayout.CENTER);
      scrollNamePanel.add(rightBtnPanel, BorderLayout.EAST);

      vis2DPanel.add(vis2DBtn, BorderLayout.CENTER);
      vis3DPanel.add(vis3DBtn, BorderLayout.CENTER);
      if(_is3D) {
	 vis2D3DPanel.add(vis2DPanel, BorderLayout.WEST);
	 vis2D3DPanel.add(vis3DPanel, BorderLayout.EAST);
      } else {
	 vis2D3DPanel.add(vis2DPanel, BorderLayout.CENTER);
      }

      colPanel.add(colBtn, BorderLayout.CENTER);
      zapPanel.add(zapBtn, BorderLayout.CENTER);
      colZapPanel.add(colPanel, BorderLayout.WEST);
      colZapPanel.add(zapPanel, BorderLayout.EAST);

      buttonPanel.add(vis2D3DPanel, BorderLayout.WEST);
      buttonPanel.add(colZapPanel, BorderLayout.EAST);

      innerPanel.add(scrollNamePanel, BorderLayout.CENTER);
      innerPanel.add(buttonPanel, BorderLayout.EAST);

      outerPanel.add(innerPanel, BorderLayout.CENTER); // set height to expand
      this.add(outerPanel, BorderLayout.CENTER); // allows text field to expand
//-----------------------------------------------------
      /* for testing
      printPanelSize(innerPanel, "outer");
      printPanelSize(outerPanel, "inner");
      printPanelSize(scrollNamePanel, "scrollName");
      printPanelSize(buttonPanel, "button");
      printPanelSize(vis2DPanel, "vis2D");
      printPanelSize(vis3DPanel, "vis3D");
      printPanelSize(vis2D3DPanel, "vis2D3D");
      printPanelSize(colPanel, "col");
      printPanelSize(zapPanel, "zap");
      printPanelSize(colZapPanel, "colZap");
      */
//......................................
      setToolTips();
//......................................

   } // makeKeyEntryGUI()

   
   private void printPanelSize(JPanel panel, String name) {
      Dimension dim = null;
      dim = panel.getPreferredSize();
      System.out.println("---------- "+name+"----------");
      System.out.println("panel W = "+dim.width);
      System.out.println("panel H = "+dim.height);
   }
//......................................

   /**
    *   Attaches tool tip text to the controls in the KeyEntry.
    */
   private void setToolTips() {

      ToolTipManager.sharedInstance().setInitialDelay( 1500 );
      ToolTipManager.sharedInstance().setDismissDelay( 1500 );

      colBtn.setToolTipText(colTipText);
      leftBtn.setToolTipText(leftTipText);
      rightBtn.setToolTipText(rightTipText);
      vis2DBtn.setToolTipText(vizTip2Text);
      vis3DBtn.setToolTipText(vizTip3Text);
      zapBtn.setToolTipText(zapTipText);

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
