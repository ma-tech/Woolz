package sectionViewer;

import java.net.*;
import java.lang.ClassLoader;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;

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

   /** The 'start' icon  */
   protected ImageIcon startIcon = null;
   /** The 'end' icon  */
   protected ImageIcon endIcon = null;
   /** The 'show' icon  */
   protected ImageIcon showIcon = null;
   /** The 'hide' icon  */
   protected ImageIcon hideIcon = null;
   /** The 'hide2D' icon  */
   protected ImageIcon hide2DIcon = null;
   /** The 'show2D' icon  */
   protected ImageIcon show2DIcon = null;
   /** The 'hide3D' icon  */
   protected ImageIcon hide3DIcon = null;
   /** The 'show3D' icon  */
   protected ImageIcon show3DIcon = null;

   /** The tool tip for the colour chooser button */
   protected String colTipText = "choose anatomy component colour";
   /** The tool tip for the left arrow button */
   protected String leftTipText = "jump to start of name";
   /** The tool tip for the right arrow button */
   protected String rightTipText = "jump to end of name";
   /** tool tip for the 2D visibility button */
   protected String show2TipText = "show 2D component";
   /** tool tip for the 3D visibility button */
   protected String show3TipText = "show 3D component";
   /** tool tip for the 2D visibility button */
   protected String hide2TipText = "hide 2D component";
   /** tool tip for the 3D visibility button */
   protected String hide3TipText = "hide 3D component";
   /** The tool tip for the visibility button */
   protected String vizTip2Text = "toggle 2D visibility";
   /** The tool tip for the visibility button */
   protected String vizTip3Text = "toggle 3D visibility";
   /** The tool tip for the zap button */
   protected String zapTipText = "remove anatomy component";
   /** The tool tip for the name Text Field */
   protected String nameTipText = "drag mouse to scroll";

   protected boolean _is3D;
   protected static int _entryW = 300;
   //protected static int _entryH = 25;
   protected static int _entryH = 20;

   protected int _buttonW1 = 17;
   protected int _buttonW2 = 25;
   protected int _buttonH = 25;
   protected int _tfW = 200;
   protected int _tfH = 25;
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

      // Get current classloader
      // this is needed to retrieve resources from a jar file (ie from web start deployed application)
      ClassLoader cl = this.getClass().getClassLoader();
      String imgPath = "images/";

      try {
        zapIcon = new ImageIcon(cl.getResource(imgPath + "bin.png"));
        startIcon = new ImageIcon(cl.getResource(imgPath + "start2.png"));
        endIcon = new ImageIcon(cl.getResource(imgPath + "end2.png"));
        showIcon = new ImageIcon(cl.getResource(imgPath + "show.png"));
        hideIcon = new ImageIcon(cl.getResource(imgPath + "hide.png"));
      }
      catch (Exception ex) {
         ex.printStackTrace();
      }

//......................................
      this.setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
      this.setPreferredSize(new Dimension(_entryW, _entryH));
      this.setMinimumSize(new Dimension(_entryW, _entryH));
      //this.setBackground(Color.green);
//......................................

      colBtn = new JButton();
      colBtn.setPreferredSize(new Dimension(_buttonW2,_entryH));
      colBtn.setMinimumSize(new Dimension(_buttonW2,_entryH));
      colBtn.setMaximumSize(new Dimension(_buttonW2,_entryH));

      textF = new ScrollableTextField("", 50);
      textF.setEditable(false);
      textF.setPreferredSize(new Dimension(_tfW,_tfH));
      textF.setMinimumSize(new Dimension(_tfW,_tfH));
      textF.setFont(new Font("default", Font.PLAIN, 10));

      leftBtn = new JButton();
      leftBtn.setPreferredSize(new Dimension(_buttonW1,_entryH));
      leftBtn.setMinimumSize(new Dimension(_buttonW1,_entryH));
      leftBtn.setMaximumSize(new Dimension(_buttonW1,_entryH));
      if(startIcon != null) {
	 leftBtn.setIcon(startIcon);
      }

      rightBtn = new JButton();
      rightBtn.setPreferredSize(new Dimension(_buttonW1,_entryH));
      rightBtn.setMinimumSize(new Dimension(_buttonW1,_entryH));
      rightBtn.setMaximumSize(new Dimension(_buttonW1,_entryH));
      if(endIcon != null) {
	 rightBtn.setIcon(endIcon);
      }

      vis2DBtn = new JButton();
      vis2DBtn.setPreferredSize(new Dimension(_buttonW2,_entryH));
      vis2DBtn.setMinimumSize(new Dimension(_buttonW2,_entryH));
      vis2DBtn.setMaximumSize(new Dimension(_buttonW2,_entryH));
      if(hideIcon != null) {
	 vis2DBtn.setIcon(hideIcon);
      }

      vis3DBtn = new JButton();
      vis3DBtn.setPreferredSize(new Dimension(_buttonW2,_entryH));
      vis3DBtn.setMinimumSize(new Dimension(_buttonW2,_entryH));
      vis3DBtn.setMaximumSize(new Dimension(_buttonW2,_entryH));
      if(hideIcon != null) {
	 vis3DBtn.setIcon(hideIcon);
      }

      zapBtn = new JButton();
      zapBtn.setPreferredSize(new Dimension(_buttonW2,_entryH));
      zapBtn.setMinimumSize(new Dimension(_buttonW2,_entryH));
      zapBtn.setMaximumSize(new Dimension(_buttonW2,_entryH));
      if(zapIcon != null) {
	 zapBtn.setIcon(zapIcon);
      }
//......................................
      this.add(leftBtn);
      this.add(rightBtn);
      this.add(textF);
      this.add(vis2DBtn);
      if(_is3D) {
	 this.add(vis3DBtn);
      }
      this.add(colBtn);
      this.add(zapBtn);
//......................................
      setToolTips();
//......................................

   } // makeKeyEntryGUI()

//......................................

   /**
    *   Attaches tool tip text to the controls in the KeyEntry.
    */
   private void setToolTips() {

      ToolTipManager.sharedInstance().setInitialDelay( 800 );
      ToolTipManager.sharedInstance().setDismissDelay( 1500 );

      colBtn.setToolTipText(colTipText);
      leftBtn.setToolTipText(leftTipText);
      rightBtn.setToolTipText(rightTipText);
      vis2DBtn.setToolTipText(hide2TipText);
      vis3DBtn.setToolTipText(hide3TipText);
      zapBtn.setToolTipText(zapTipText);
      textF.setToolTipText(nameTipText);

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
