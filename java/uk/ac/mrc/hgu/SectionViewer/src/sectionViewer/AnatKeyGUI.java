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
 *   The GUI for AnatKey.
 */
public class AnatKeyGUI extends JFrame {

   /**
    *   The system file separator ("/" or "\").
    */
   private String SLASH = System.getProperty("file.separator");

   /**
    *   The number of distinct colours for the AnatKey.
    */
   protected final int _ncols = 6;

   /**
    *   The initial set of colours for AnatKey entries.
    */
   public static int _cols[];

   /**
    *   The set of red, green, blue and transparency values
    *   which generate the initial colours for the AnatKey entries.
    */
   protected int _rgbt[][] = {{128,96,0,125},
                              {255,0,0,125},
                              {255,255,0,125},
                              {0,0,255,125},
                              {0,255,0,125},
                              {0,255,255,125}};

   /**
    *   The 'Content Pane' of this component
    */
   protected JPanel topPanel = null;

   /**
    *   The container for rows of the key.
    */
   protected JPanel kTopPanel = null;

   /**
    *   The row header.
    */
   protected KeyHeader _header = null;

   /**
    *   The initial width of the anatomy key.
    */
   protected int _keyWidth = 400;

   /**
    *   The empty height of the anatomy key.
    */
   protected final int _emptyH = 75;

   protected static boolean _is3D;
//-------------------------------------------------------------
   /**
    *   Constructor. 
    */
   protected AnatKeyGUI(String str) {
      this(str, false);
   }
   
   /**
    *   Constructor. 
    */
   protected AnatKeyGUI(String str, boolean is3D) {
      super(str);
      _is3D = is3D;
      setCols();
      makeGUI();
   }
   
//-------------------------------------------------------------
   /**
    *   Initialises the AnatKey with its default colours.
    */
   protected void setCols() {

      _cols = new int[_ncols];
      _cols[0] = (_rgbt[0][2])|
                 (_rgbt[0][1]<<8)|
                 (_rgbt[0][0]<<16)|
                 (_rgbt[0][3]<<24);
      _cols[1] = (_rgbt[1][2])|
                 (_rgbt[1][1]<<8)|
                 (_rgbt[1][0]<<16)|
                 (_rgbt[1][3]<<24);
      _cols[2] = (_rgbt[2][2])|
                 (_rgbt[2][1]<<8)|
                 (_rgbt[2][0]<<16)|
                 (_rgbt[2][3]<<24);
      _cols[3] = (_rgbt[3][2])|
                 (_rgbt[3][1]<<8)|
                 (_rgbt[3][0]<<16)|
                 (_rgbt[3][3]<<24);
      _cols[4] = (_rgbt[4][2])|
                 (_rgbt[4][1]<<8)|
                 (_rgbt[4][0]<<16)|
                 (_rgbt[4][3]<<24);
      _cols[5] = (_rgbt[5][2])|
                 (_rgbt[5][1]<<8)|
                 (_rgbt[5][0]<<16)|
                 (_rgbt[5][3]<<24);

   }
//-------------------------------------------------------------
   /**
    *   Called by the constructor. Most of the work of building the GUI
    *   is done here.
    */
   private void makeGUI() {

      topPanel = (JPanel)this.getContentPane();
      topPanel.setLayout(new BorderLayout());
      //topPanel.setBackground(Color.yellow);

      _header = new KeyHeader(_is3D);

      kTopPanel = new JPanel();
      kTopPanel.setLayout(new BoxLayout(kTopPanel, BoxLayout.Y_AXIS));
      kTopPanel.setPreferredSize(new Dimension(_keyWidth,0));
      //kTopPanel.setBackground(Color.red);

      topPanel.add(_header, BorderLayout.NORTH);
      topPanel.add(kTopPanel, BorderLayout.CENTER);
      
//......................................
   }

} // class AnatKeyGUI
