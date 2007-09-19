package sectionViewer;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 *   Header for AnatomyKey.
 */
public class KeyHeader extends JPanel{


   /** The system file separator ("/" or "\"). */
   private String SLASH = System.getProperty("file.separator");

   protected JLabel colLabel = null;
   protected JLabel vizLabel = null;
   protected JLabel viz2DLabel = null;
   protected JLabel viz3DLabel = null;
   protected JLabel zapLabel = null;
   protected JLabel nameLabel = null;

   protected final Font _labelFont1 = new Font("Default",Font.PLAIN,14);
   protected final Font _labelFont2 = new Font("Default",Font.PLAIN,9);
   protected final Font _labelFont3 = new Font("Default",Font.PLAIN,8);

   protected double _rads = 1.5705;
   protected boolean _is3D;

   protected int _headerW = 600;
   protected int _headerH = 20;
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

      this.setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
      this.setBackground(new Color(200, 200, 205));
      this.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
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

      nameLabel = new JLabel("Anatomy Name");
      nameLabel.setFont(_labelFont1);
      nameLabel.setHorizontalAlignment(JLabel.CENTER);
      nameLabel.setVerticalAlignment(JLabel.CENTER);

      this.add(Box.createRigidArea(new Dimension(50,0)));
      this.add(Box.createHorizontalGlue());
      this.add(nameLabel);
      this.add(Box.createHorizontalGlue());
      this.add(viz2DLabel);
      if(_is3D) {
	 this.add(Box.createRigidArea(new Dimension(13,0)));
	 this.add(viz3DLabel);
      }
      this.add(Box.createRigidArea(new Dimension(13,0)));
      this.add(colLabel);
      this.add(Box.createRigidArea(new Dimension(10,0)));
      this.add(zapLabel);
      //this.add(Box.createRigidArea(new Dimension(3,0)));
      this.add(Box.createRigidArea(new Dimension(20,0)));

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
   public int getH() {
      return _headerH;
   }

   /**
    *   Sets the width of the header
    *   @param w the width to be set.
    */
   public void setW(int w) {
      _headerW = w;
   }

} // class KeyHeader
