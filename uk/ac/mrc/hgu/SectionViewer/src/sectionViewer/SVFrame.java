package sectionViewer;

import java.awt.*;
import javax.swing.*;

/**
 *   Associates a SectionViewer with JFrame.
 */
public class SVFrame extends JFrame {

   private final boolean _debug = false;

   /**   the SectionViewer that is associated with this SVFrame */
   protected SectionViewer _SV = null;
   JPanel _contentPane;
   JRootPane _rootPane;

   /**   the initial height of the window. */
   int totalH = 500;
   /**   the initial width of the window. */
   int totalW = 400;
   int pad = 1;

   //=========================================================
   // constructor
   //=========================================================
   /**
    *   Constructs an SVFrame with the given title.
    *   @param viewstr
    */
   public SVFrame(String viewstr) {

      super(viewstr);

      if (_debug) System.out.println("enter SVFrame");

      try {
	 guiInit();
      }
      catch(Exception e) {
	 System.out.println("SVFrame");
	 e.printStackTrace();
      }

      if (_debug)
	 System.out.println("exit SVFrame");
   } // constructor

   //---------------------------------------------
   /**
    *   Builds the SVFrame.
    */
   private void guiInit() throws Exception {

      _contentPane = (JPanel)this.getContentPane();
      _rootPane = this.getRootPane();

      _contentPane.setLayout(new BorderLayout());
      _contentPane.setPreferredSize(new Dimension(totalW+pad, totalH+pad));
      _contentPane.setMinimumSize(new Dimension(totalW+pad, totalH+pad));
   }

   //---------------------------------------------
   /**
    *   Returns the SectionViewer for this SVFrame.
    */
   public SectionViewer getSectionViewer() {
      return _SV;
   }

   //---------------------------------------------
   /**
    *   Sets the SectionViewer for this SVFrame.
    *   @param sv the SectionViewer to be associated with this SVFrame.
    */
   public void setSectionViewer(SectionViewer sv) {
      _SV = sv;
      _contentPane.add(_SV, BorderLayout.CENTER);
   }

} // class SVFrame
