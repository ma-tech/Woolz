package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import javax.swing.*;
import javax.swing.border.*;
import java.io.*;

public class SVFrame extends JFrame {

   private final boolean _debug = false;

   protected SectionViewer _SV = null;
   JPanel _contentPane;
   JRootPane _rootPane;

   int totalH = 500;
   int totalW = 400;
   int pad = 1;

   //=========================================================
   // constructor
   //=========================================================
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
   private void guiInit() throws Exception {

      _contentPane = (JPanel)this.getContentPane();
      _rootPane = this.getRootPane();

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

} // class SVFrame
