package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;
import java.io.*;

public class SVFrameInt extends JInternalFrame {

   private final boolean _debug = false;

   protected SectionViewer _SV = null;
   JPanel _contentPane;

   int totalH = 500;
   int totalW = 400;
   int pad = 1;

   //=========================================================
   // constructor
   //=========================================================
   public SVFrameInt(String viewstr) {

      super(viewstr);

      if (_debug) System.out.println("enter SVFrameInt");

      try {
	 guiInit();
      }
      catch(Exception e) {
	 System.out.println("SVFrameInt");
	 e.printStackTrace();
      }

      if (_debug)
	 System.out.println("exit SVFrameInt");
   } // constructor

   //---------------------------------------------
   private void guiInit() throws Exception {

      _contentPane = (JPanel)this.getContentPane();

      _contentPane.setLayout(new BorderLayout());
      _contentPane.setPreferredSize(new Dimension(totalW+pad, totalH+pad));
      _contentPane.setMinimumSize(new Dimension(totalW+pad, totalH+pad));
   }

   //---------------------------------------------
   public SectionViewer getSectionViewer() {
      return _SV;
   }

   //---------------------------------------------
   public void setSectionViewer(SectionViewer sv) {
      _SV = sv;
      _contentPane.add(_SV, BorderLayout.CENTER);
   }

} // class SVExt
