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
   protected JLabel scrollLabel = null;
   protected JLabel vizLabel = null;
   protected JLabel viz2DLabel = null;
   protected JLabel viz3DLabel = null;
   protected JLabel zapLabel = null;
   protected JLabel nameLabel = null;

   protected String _fontFamilies[] = null;
   protected Font _fonts[] = null;

   protected Font _defFont = null;
   protected Font _labelFont1 = null;
   protected Font _labelFont2 = null;
   protected Font _labelFont3 = null;

   protected double _rads = 1.5705;
   protected boolean _is3D;

   protected static final int _headerW = 200;
   protected static final int _headerH = 40;
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
      setFonts();
      makeGUI();
   }
   
//-------------------------------------------------------------
   /**
    *   Set the font(s) for the header.
    */
   private void setFonts() {

      getFontNames(false);

      _defFont = getFont();
      /*
      System.out.println("default font name "+_defFont.getName());
      System.out.println("size in points "+
          Float.toString(_defFont.getSize2D()));
      System.out.println("plain "+_defFont.isPlain());
      System.out.println("bold "+_defFont.isBold());
      System.out.println("italic "+_defFont.isItalic());
      */

      _labelFont1 = new Font("Dialog",Font.PLAIN,16);
      _labelFont2 = new Font("Dialog",Font.PLAIN,9);
      _labelFont3 = new Font("Dialog",Font.BOLD,10);

   }
//-------------------------------------------------------------
   /**
    *   Get the font(s) available to Java.
    */
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
//-------------------------------------------------------------
   /**
    *   Called by the constructor. All of the work of building
    *   the header is done here.
    */
   private void makeGUI() {

      boolean hasAlpha = false;

      int W1 = 20;
      int W2 = 25;
      int W3 = 30;
      int H1 = 15;
      int H2 = 20;
      int H3 = 35;
      int gap = 1;

//......................................
      this.setLayout(new BorderLayout());
      this.setPreferredSize(new Dimension(0, _headerH));
      //this.setBackground(Color.green);
      this.setBorder(BorderFactory.createEtchedBorder(
                                        EtchedBorder.RAISED));
//......................................
      colLabel = new JLabel("Change");
      colLabel.setFont(_labelFont2);
      colLabel.setUI( new VLabelUI(false) );
      colLabel.setHorizontalAlignment(JLabel.CENTER);
      colLabel.setVerticalAlignment(JLabel.CENTER);
      scrollLabel = new JLabel("Scroll");
      scrollLabel.setFont(_labelFont3);
      /*
      scrollLabel.setHorizontalAlignment(JLabel.CENTER);
      scrollLabel.setVerticalAlignment(JLabel.CENTER);
      */
      vizLabel = new JLabel("Visibility");
      vizLabel.setFont(_labelFont3);
      vizLabel.setHorizontalAlignment(JLabel.CENTER);
      viz2DLabel = new JLabel("2D");
      viz2DLabel.setFont(_labelFont3);
      viz2DLabel.setHorizontalAlignment(JLabel.CENTER);
      viz3DLabel = new JLabel("3D");
      viz3DLabel.setFont(_labelFont3);
      viz3DLabel.setHorizontalAlignment(JLabel.CENTER);
      //zapLabel = new JLabel("Delete");
      zapLabel = new JLabel("Remove");
      zapLabel.setFont(_labelFont2);
      zapLabel.setUI( new VLabelUI(false) );
      zapLabel.setHorizontalAlignment(JLabel.CENTER);
      zapLabel.setVerticalAlignment(JLabel.CENTER);
      nameLabel = new JLabel("Full Component Name");
      nameLabel.setFont(_labelFont1);
      nameLabel.setHorizontalAlignment(JLabel.CENTER);
      nameLabel.setVerticalAlignment(JLabel.CENTER);
//......................................
      /**
       *   Outer container for header
       *   which maintains the preferred size.
       *   <br>Contains kPanel0.
       */
      JPanel kPanel = new JPanel();
      kPanel.setLayout(new BorderLayout(gap,gap));
      //kPanel.setBackground(Color.orange);
      kPanel.setPreferredSize(new Dimension(_headerW,_headerH));

      /**
       *   Inner container for header.
       *   <br>Contains kPanel00, txtPanel0, kPanel02.
       */
      JPanel kPanel0 = new JPanel();
      kPanel0.setLayout(new BorderLayout(gap,gap));
      //kPanel0.setBackground(Color.pink);
      kPanel0.setPreferredSize(new Dimension(0,_headerH));

      /**
       *   Outer container for <em>colour chooser</em>
       *   and <em>scroll</em> labels.
       *   <br>Contains colPanel0, navPanel0.
       *   <br>Added to kPanel0.
       */
      JPanel kPanel00 = new JPanel();
      kPanel00.setLayout(new BorderLayout(gap,gap));
      kPanel00.setPreferredSize(new Dimension(2*W2,0));
      kPanel00.setBackground(Color.blue);
      /*
      kPanel00.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */

      /**
       *   Outer container for <em>visibility</em> labels
       *   <br>Contains
       *   <br>Added to kPanel02.
       */
      JPanel kPanel01 = new JPanel();
      kPanel01.setLayout(new BorderLayout(gap,gap));
      kPanel01.setPreferredSize(new Dimension(2*W2,0));
      //kPanel01.setBackground(Color.magenta);
      kPanel01.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      /**
       *   Outer container for <em>visibility</em> and
       *   <em>zap</em> labels.
       *   <br>Contains kPanel02, btnPanel03.
       *   <br>Added to kPanel0.
       */
      JPanel kPanel02 = new JPanel();
      kPanel02.setLayout(new BorderLayout(gap,gap));
      kPanel02.setPreferredSize(new Dimension(W1+2*W2,0));
      kPanel02.setBackground(Color.red);
      /*
      kPanel02.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */

      /**
       *   Container for <em>visibility</em> label.
       *   <br>Added to kPanel01.
       */
      JPanel kPanel03 = new JPanel();
      kPanel03.setLayout(new BorderLayout(gap,gap));
      kPanel03.setPreferredSize(new Dimension(0,H1));
      //kPanel02.setBackground(Color.red);
      /*
      kPanel02.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */

      /**
       *   Container for <em>2D</em> and <em>3D</em> labels.
       *   <br>Added to kPanel01.
       */
      JPanel kPanel04 = new JPanel();
      kPanel04.setLayout(new BorderLayout(gap,gap));
      kPanel04.setPreferredSize(new Dimension(0,H2));
      //kPanel02.setBackground(Color.red);
      /*
      kPanel02.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */

      JPanel colPanel0 = new JPanel();
      colPanel0.setLayout(new BorderLayout(gap,gap));
      colPanel0.setPreferredSize(new Dimension(W1-4*gap,0));
      colPanel0.add(colLabel);

      JPanel txtPanel0 = new JPanel();
      txtPanel0.setLayout(new BorderLayout(gap,gap));
      txtPanel0.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      txtPanel0.add(nameLabel, BorderLayout.CENTER);

      JPanel navPanel0 = new JPanel();
      navPanel0.setLayout(new BorderLayout(gap,gap));
      navPanel0.setPreferredSize(new Dimension(W3+4*gap,0));
      navPanel0.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      navPanel0.add(scrollLabel, BorderLayout.CENTER);

      JPanel btnPanel02 = new JPanel();
      btnPanel02.setPreferredSize(new Dimension(W2,0));
      btnPanel02.setLayout(new BorderLayout(gap,gap));
      /*
      btnPanel02.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */
      btnPanel02.add(viz2DLabel, BorderLayout.CENTER);

      JPanel btnPanel03 = new JPanel();
      btnPanel03.setPreferredSize(new Dimension(W2,0));
      btnPanel03.setLayout(new BorderLayout(gap,gap));
      /*
      btnPanel03.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */
      btnPanel03.add(viz3DLabel, BorderLayout.CENTER);

      JPanel btnPanel04 = new JPanel();
      btnPanel04.setPreferredSize(new Dimension(W1,0));
      btnPanel04.setLayout(new BorderLayout(gap,gap));
      /*
      btnPanel04.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      */
      btnPanel04.add(zapLabel, BorderLayout.CENTER);

      kPanel00.add(colPanel0, BorderLayout.WEST);
      kPanel00.add(navPanel0, BorderLayout.EAST);
      kPanel01.add(kPanel03, BorderLayout.NORTH);
      kPanel01.add(kPanel04, BorderLayout.SOUTH);
      kPanel02.add(kPanel01, BorderLayout.WEST);
      kPanel02.add(btnPanel04, BorderLayout.EAST);
      kPanel03.add(vizLabel, BorderLayout.CENTER);
      if(_is3D) {
	 kPanel04.add(btnPanel02, BorderLayout.WEST);
	 kPanel04.add(btnPanel03, BorderLayout.EAST);
      } else {
	 kPanel04.add(btnPanel02, BorderLayout.CENTER);
      }
      kPanel0.add(kPanel00, BorderLayout.WEST);
      kPanel0.add(txtPanel0, BorderLayout.CENTER);
      kPanel0.add(kPanel02, BorderLayout.EAST);

      //kPanel.add(kPanel0, BorderLayout.NORTH); // set height to preferred size
      kPanel.add(kPanel0, BorderLayout.CENTER); // set height to expand
      this.add(kPanel, BorderLayout.CENTER); // allows text field to expand

   } // makeGUI()

//-------------------------------------------------------------
   /**
    *   Returns the height of the header
    *   @return _headerH.
    */
   public static int getH() {
      return _headerH;
   }

} // class KeyHeader
