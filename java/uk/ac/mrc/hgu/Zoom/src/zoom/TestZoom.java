import zoom.*;

import javax.swing.*;
import java.awt.*;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;


public class TestZoom extends JFrame {

   static JFrame frame;
   JPanel topPanel = null;

   ZoomToZoomAdaptor ZZ = null;

   public TestZoom() {
      super("TestZoom");

      topPanel = (JPanel)this.getContentPane();
      topPanel.setLayout(new BorderLayout());
      Zoom zoom = new Zoom();
      ZZ = new ZoomToZoomAdaptor(zoom);
      zoom.addChangeListener(ZZ);
      topPanel.add(zoom, BorderLayout.NORTH);
      zoom.setEnabled(true);
   }

//---------------------------------------
  // change 'zoom' using Zoom control
  public class ZoomToZoomAdaptor
      implements ChangeListener {

    Zoom control;
    int val;

    public ZoomToZoomAdaptor(Zoom cntrl) {
      control = cntrl;
    }

    public void stateChanged(ChangeEvent e) {
      val = control.getValue();
      System.out.println("ZoomToZoomAdaptor: stateChanged "+val);
    }
  }


  public static void main(String[] args) {
    frame = new TestZoom();
    frame.setTitle("Test");
    frame.setSize(200, 100);
    frame.setVisible(true);
    frame.addWindowListener(new WindowAdapter() {
      public void windowClosing(WindowEvent e) {
        System.exit(0);
      }
    });
  }

}
