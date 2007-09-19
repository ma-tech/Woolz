package sectionViewer;

import java.awt.*;
import javax.swing.*;

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
   //protected final int _ncols = 6;
   protected final int _ncols = 12;

   /**
    *   The initial set of colours for AnatKey entries.
    */
   public static int _cols[];

   /**
    *   The set of red, green, blue and transparency values
    *   which generate the initial colours for the AnatKey entries.
    */
   protected int _rgbt[][] = {
                              {255,0,0,125},     // red
                              {0,255,0,125},     // green
                              {0,0,255,125},     // blue
                              {255,255,0,125},   // yellow
                              {0,255,255,125},   // cyan
                              {255,0,255,125},   // magenta
                              {128,96,0,125},    // brown
			      {0,108,0,125},        // dark green
			      {100,149,237,125},   // cornflower blue
                              {255,96,0,125},    // orange
                              {190,220,220,125}, // grey
                              {255,222,173,125}    // navajo white
			      };

   /**
    *   The 'Content Pane' of this component
    */
   protected JPanel topPanel = null;

   /**
    *   The container for rows of the key.
    */
   protected ScrollablePanel kTopPanel = null;

   /**
    *   ScrollPane for kTopPanel
    */
   protected JScrollPane kScroll = null;

   /**
    *   The row header.
    */
   protected KeyHeader _header = null;

   protected static boolean _is3D;
   //protected int _goldenW = 600;
   //protected int _goldenH = 371;
   protected int _goldenW = 500;
   protected int _goldenH = 309;
   private int _minW = 350;
   private int _minH = 30;

   //protected Image keyImage = null;
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
      _cols[6] = (_rgbt[6][2])|
                 (_rgbt[6][1]<<8)|
                 (_rgbt[6][0]<<16)|
                 (_rgbt[6][3]<<24);
      _cols[7] = (_rgbt[7][2])|
                 (_rgbt[7][1]<<8)|
                 (_rgbt[7][0]<<16)|
                 (_rgbt[7][3]<<24);
      _cols[8] = (_rgbt[8][2])|
                 (_rgbt[8][1]<<8)|
                 (_rgbt[8][0]<<16)|
                 (_rgbt[8][3]<<24);
      _cols[9] = (_rgbt[9][2])|
                 (_rgbt[9][1]<<8)|
                 (_rgbt[9][0]<<16)|
                 (_rgbt[9][3]<<24);
      _cols[10] = (_rgbt[10][2])|
                 (_rgbt[10][1]<<8)|
                 (_rgbt[10][0]<<16)|
                 (_rgbt[10][3]<<24);
      _cols[11] = (_rgbt[11][2])|
                 (_rgbt[11][1]<<8)|
                 (_rgbt[11][0]<<16)|
                 (_rgbt[11][3]<<24);

   }
//-------------------------------------------------------------
   /**
    *   Called by the constructor. Most of the work of building the GUI
    *   is done here.
    */
   private void makeGUI() {

      // Get current classloader
      // this is needed to retrieve resources from a jar file (ie from web start deployed application)
      ClassLoader cl = this.getClass().getClassLoader();
      String imgPath = "images/";
      ImageIcon imgIcon = new ImageIcon(cl.getResource(imgPath + "AV.png"));
      if (imgIcon != null) {
	 setIconImage(imgIcon.getImage());
      }

      this.setResizable(true);
      this.setPreferredSize(new Dimension(_goldenW, _goldenH));

      topPanel = (JPanel)this.getContentPane();
      topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.PAGE_AXIS));

      _header = new KeyHeader(_is3D);

      kTopPanel = new ScrollablePanel();
      kTopPanel.setLayout(new BoxLayout(kTopPanel, BoxLayout.PAGE_AXIS));

      kScroll = new JScrollPane(kTopPanel);
      kScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
      kScroll.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      kScroll.setMinimumSize(new Dimension(_minW, _minH));

      topPanel.add(_header);
      topPanel.add(kScroll);

      this.pack();

   }
//......................................
   class ScrollablePanel extends JPanel implements Scrollable {
      
      public ScrollablePanel() {
         super();
      }
      public Dimension getPreferredScrollableViewportSize() {
         int height = this.getPreferredSize().height;
         int width = this.getPreferredSize().width;
	 
	 return new Dimension(width, height);
      }
      public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
         return 1;
      }
      public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
         return 1;
      }
      public boolean getScrollableTracksViewportWidth() {
         return true;
      }
      public boolean getScrollableTracksViewportHeight() {
         return false;
      }
   }

} // class AnatKeyGUI
