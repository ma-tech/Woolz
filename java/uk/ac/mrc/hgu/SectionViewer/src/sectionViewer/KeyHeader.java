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
 *   Header for AnatomyKey.
 */
public class KeyHeader extends JPanel{


   /** The system file separator ("/" or "\"). */
   private String SLASH = System.getProperty("file.separator");

   protected JLabel colLabel = null;
   //protected JLabel scrollLabel = null;
   protected JLabel vizLabel = null;
   protected JLabel viz2DLabel = null;
   protected JLabel viz3DLabel = null;
   protected JLabel zapLabel = null;
   protected JLabel nameLabel = null;

   //protected String _fontFamilies[] = null;
   //protected Font _fonts[] = null;
   //protected Font _defFont = null;

   protected final Font _labelFont1 = new Font("Default",Font.PLAIN,14);
   protected final Font _labelFont2 = new Font("Default",Font.PLAIN,9);
   protected final Font _labelFont3 = new Font("Default",Font.PLAIN,8);

   protected double _rads = 1.5705;
   protected boolean _is3D;

   protected static int _headerW = 200;
   protected static int _headerH;
//-------------------------------------------------------------
   /**
    *   Constructor. 
    */
   protected KeyHeader() {
      this(false);
   }
   
   /**
    *   Constructor. 
    */
   protected KeyHeader(boolean is3D) {
      _is3D = is3D;
      if(_is3D) {
         _headerH = 32;
      } else {
         _headerH = 30;
      }
      makeGUI();
   }
   
//-------------------------------------------------------------
   /*
   private void getFontNames(boolean bool) {
      _fontFamilies = GraphicsEnvironment.getLocalGraphicsEnvironment().
                 getAvailableFontFamilyNames();

      _fonts = GraphicsEnvironment.getLocalGraphicsEnvironment().
                 getAllFonts();

      if(bool) {
	 for(int i=0; i<_fontFamilies.length; i++) {
	    System.out.println(_fonts[i]);
	 }
	 for(int j=0; j<_fonts.length; j++) {
	    System.out.println(_fonts[j].getFontName());
	 }
      }
   }
   */
//-------------------------------------------------------------
   /**
    *   Called by the constructor. All of the work of building
    *   the header is done here.
    */
   private void makeGUI() {

      boolean hasAlpha = false;

      int btnW1 = 18; // left, right button panels
      int btnW2 = 20; // colour, viz & zap button panels
      int btnW3 = 25; // combined left & right button panel

      int W1 = 20;
      int W2 = 25;
      int W3 = 30;
      int H1 = 10;
      int H2 = 20;
      int H3 = 35;
      int gap = 1;

      int w = 0;

//......................................
      this.setLayout(new BorderLayout());
      this.setPreferredSize(new Dimension(0, _headerH));
      //this.setBackground(Color.green);
      this.setBorder(BorderFactory.createEtchedBorder(
                                        EtchedBorder.RAISED));
//......................................
      colLabel = new JLabel("Col");
      colLabel.setFont(_labelFont2);
      colLabel.setHorizontalAlignment(JLabel.CENTER);
      colLabel.setVerticalAlignment(JLabel.CENTER);
      vizLabel = new JLabel("See");
      vizLabel.setFont(_labelFont2);
      vizLabel.setHorizontalAlignment(JLabel.CENTER);
      viz2DLabel = new JLabel("2D");
      viz2DLabel.setFont(_labelFont3);
      viz2DLabel.setHorizontalAlignment(JLabel.CENTER);
      viz2DLabel.setVerticalAlignment(JLabel.BOTTOM);
      viz3DLabel = new JLabel("3D");
      viz3DLabel.setFont(_labelFont3);
      viz3DLabel.setHorizontalAlignment(JLabel.CENTER);
      viz3DLabel.setVerticalAlignment(JLabel.BOTTOM);
      zapLabel = new JLabel("Del");
      zapLabel.setFont(_labelFont2);
      zapLabel.setHorizontalAlignment(JLabel.CENTER);
      zapLabel.setVerticalAlignment(JLabel.CENTER);
      nameLabel = new JLabel("Anatomy Component");
      nameLabel.setFont(_labelFont1);
      nameLabel.setHorizontalAlignment(JLabel.CENTER);
      nameLabel.setVerticalAlignment(JLabel.CENTER);
//KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK

      /**
       *   Outer container for header
       *   which maintains the preferred size.
       *   <br>Contains innerPanel.
       */
      JPanel outerPanel = new JPanel();
      outerPanel.setLayout(new BorderLayout(gap,gap));
      //outerPanel.setBackground(Color.orange);
      outerPanel.setPreferredSize(new Dimension(_headerW,_headerH));

      /**
       *   Inner container for header.
       */
      JPanel innerPanel = new JPanel();
      innerPanel.setLayout(new BorderLayout(gap,gap));
      //innerPanel.setBackground(Color.pink);
      innerPanel.setPreferredSize(new Dimension(0,_headerH));

//....................................
      /**
       *   Container for <em>component name</em> label
       *   <br>Added to innerPanel.
       */
      JPanel scrollNamePanel = new JPanel();
      scrollNamePanel.setLayout(new BorderLayout(gap,gap));
      scrollNamePanel.setPreferredSize(new Dimension(0,_headerH));
      scrollNamePanel.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
//....................................

      /**
       *   Container for <em>visibility</em> label.
       *   <br>Added to visPanel.
       */
      JPanel visVisPanel = new JPanel();
      visVisPanel.setLayout(new BorderLayout(gap,gap));

//....................................
      /**
       *   Container for <em>vis2D</em> label.
       *   <br>Added to vis2D3DBtnPanel.
       *   Only used for 3D case.
       */
      JPanel vis2DBtnPanel = new JPanel();
      vis2DBtnPanel.setLayout(new BorderLayout(gap,gap));
      vis2DBtnPanel.setPreferredSize(new Dimension(btnW3, H2));

      /**
       *   Container for <em>vis3D</em> label.
       *   <br>Added to vis2D3DBtnPanel.
       *   Only used for 3D case.
       */
      JPanel vis3DBtnPanel = new JPanel();
      vis3DBtnPanel.setLayout(new BorderLayout(gap,gap));
      vis3DBtnPanel.setPreferredSize(new Dimension(btnW3, H2));

      /**
       *   Container for <em>2D,3D</em> labels.
       *   <br>Added to visPanel.
       *   Only used for 3D case.
       */
      JPanel vis2D3DBtnPanel = new JPanel();
      vis2D3DBtnPanel.setLayout(new BorderLayout(gap,gap));
      vis2D3DBtnPanel.setPreferredSize(new Dimension(2*btnW3, H2));
//....................................

      /**
       *   Overall container for <em>visibility</em> labels
       *   <br>Added to buttonPanel.
       */
      JPanel visPanel = new JPanel();
      visPanel.setLayout(new BorderLayout(gap,gap));
      visPanel.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      //visPanel.setBackground(Color.magenta);

//....................................

      JPanel colBtnPanel = new JPanel();
      colBtnPanel.setLayout(new BorderLayout(gap,gap));
      colBtnPanel.setPreferredSize(new Dimension(btnW3+2*gap, _headerH));
      colBtnPanel.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      //colBtnPanel.setBackground(Color.pink);

      JPanel zapBtnPanel = new JPanel();
      zapBtnPanel.setPreferredSize(new Dimension(btnW3-2*gap, _headerH));
      zapBtnPanel.setLayout(new BorderLayout(gap,gap));
      zapBtnPanel.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      //zapBtnPanel.setBackground(Color.yellow);

      JPanel colZapBtnPanel = new JPanel();
      w = colBtnPanel.getPreferredSize().width;
      w += zapBtnPanel.getPreferredSize().width;
      colZapBtnPanel.setPreferredSize(new Dimension(w, _headerH));
      w = 0;
      colZapBtnPanel.setLayout(new BorderLayout(gap,gap));
      //colZapBtnPanel.setBackground(Color.green);

//....................................
      /**
       *   Container for <em>visibility</em>
       *   <em>col</em> and <em>zap</em> labels.
       *   <br>Added to innerPanel.
       */
      JPanel buttonPanel = new JPanel();
      buttonPanel.setLayout(new BorderLayout(gap,gap));
      buttonPanel.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      //buttonPanel.setBackground(Color.red);

//KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      scrollNamePanel.add(nameLabel, BorderLayout.CENTER);

      visVisPanel.add(vizLabel, BorderLayout.CENTER);
      vis2DBtnPanel.add(viz2DLabel, BorderLayout.SOUTH);
      vis3DBtnPanel.add(viz3DLabel, BorderLayout.SOUTH);
      vis2D3DBtnPanel.add(vis2DBtnPanel, BorderLayout.WEST);
      vis2D3DBtnPanel.add(vis3DBtnPanel, BorderLayout.EAST);

      colBtnPanel.add(colLabel, BorderLayout.CENTER);
      zapBtnPanel.add(zapLabel, BorderLayout.CENTER);
      colZapBtnPanel.add(colBtnPanel, BorderLayout.WEST);
      colZapBtnPanel.add(zapBtnPanel, BorderLayout.EAST);

      if(_is3D) {
         visVisPanel.setPreferredSize(new Dimension(2*btnW3, H1));
         visPanel.setPreferredSize(new Dimension(2*btnW3, _headerH));
	 visPanel.add(visVisPanel, BorderLayout.NORTH);
	 visPanel.add(vis2D3DBtnPanel, BorderLayout.SOUTH);
      } else {
         visVisPanel.setPreferredSize(new Dimension(btnW3, _headerH));
         visPanel.setPreferredSize(new Dimension(btnW3+2*gap, _headerH));
	 visPanel.add(visVisPanel, BorderLayout.CENTER);
      }
      w = colZapBtnPanel.getPreferredSize().width;
      w += visPanel.getPreferredSize().width;
      buttonPanel.setPreferredSize(new Dimension(w, 0));
      w = 0;

      buttonPanel.add(visPanel, BorderLayout.WEST);
      buttonPanel.add(colZapBtnPanel, BorderLayout.EAST);

      innerPanel.add(scrollNamePanel, BorderLayout.CENTER);
      innerPanel.add(buttonPanel, BorderLayout.EAST);

      outerPanel.add(innerPanel, BorderLayout.CENTER); // set height to expand
      this.add(outerPanel, BorderLayout.CENTER); // allows text field to expand
      //KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
      /* for testing
      printPanelSize(innerPanel, "outer");
      printPanelSize(outerPanel, "inner");
      printPanelSize(scrollNamePanel, "scrollName");
      printPanelSize(buttonPanel, "button");
      printPanelSize(visVisPanel, "visVis");
      printPanelSize(vis2DBtnPanel, "vis2DBtn");
      printPanelSize(vis3DBtnPanel, "vis3DBtn");
      printPanelSize(visPanel, "vis");
      printPanelSize(vis2D3DBtnPanel, "vis2D3DBtn");
      printPanelSize(colBtnPanel, "colBtn");
      printPanelSize(zapBtnPanel, "zapBtn");
      printPanelSize(colZapBtnPanel, "colZapBtn");
      */

   } // makeGUI()

   private void printPanelSize(JPanel panel, String name) {
      Dimension dim = null;
      dim = panel.getPreferredSize();
      System.out.println("---------- "+name+"----------");
      System.out.println("panel W = "+dim.width);
      System.out.println("panel H = "+dim.height);
   }
//-------------------------------------------------------------
   /**
    *   Returns the height of the header
    *   @return _headerH.
    */
   public static int getH() {
      return _headerH;
   }

} // class KeyHeader
