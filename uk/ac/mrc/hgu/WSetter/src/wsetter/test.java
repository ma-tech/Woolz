/*!
* \file         test.java
* \author       Nick Burton
* \date         April 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        test prog for WSetter
* \todo         -
* \bug          None known.
*/
package wsetter;

import wsetter.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

public class test extends JFrame{

static private String ident = "MRC HGU $Id$";

// VIEW ADAPTOR
//---------------------------------------
  public static class ModelToPanelAdaptor
    implements ChangeListener {
    WSetter model;
    JPanel view;

    public ModelToPanelAdaptor(WSetter mdl, JPanel vw) {
      model = mdl;
      view = vw;
    }

    public void stateChanged(ChangeEvent e) {
      int col = (int)model.getValue();
      view.setBackground(new Color(col,100,200));
    }
  } // class sliderToModelAdaptor
//---------------------------------------

  public test(String title) {
   super(title);
   // top panel to hold the bits
   JPanel topPanel = new JPanel();
   topPanel.setLayout(new FlowLayout());
   topPanel.setBackground(new Color(255,0,255));

   // the combi-sliders
   WSetter yawSetter1 = new WSetter();
   yawSetter1.setSliderLabel("Yaw");
   yawSetter1.setBgc(new Color(240, 250, 240));
   yawSetter1.setModelMax(250);
   yawSetter1.setSlidingEvents(false);
   yawSetter1.setModelInitVal(20.0f);
   yawSetter1.setLabelWidth(50);
   System.out.println("value = "+ yawSetter1.getValue() + "\n");

   WSetter yawSetter2 = new WSetter();
   yawSetter2.setSliderLabel("Roll");
   yawSetter2.setBgc(new Color(240, 240, 250));
   yawSetter2.setModelMin(10);
   yawSetter2.setModelMax(210);
   yawSetter2.setSlidingEvents(true);
   yawSetter2.setWidth(500);
   yawSetter2.setHeight(30);
   yawSetter2.setModelInitVal(200.0f);
   yawSetter2.setTextWidth(80);
   yawSetter2.setValue(50.0f);

   // something to change
   JPanel viewPanel = new JPanel();
   viewPanel.setPreferredSize(new Dimension(100,100));

   // add the bits
   topPanel.add(yawSetter1);
   topPanel.add(yawSetter2);
   topPanel.add(viewPanel);

   this.getContentPane().add(topPanel);

    // adapters for MVC model
    ModelToPanelAdaptor MToP = new ModelToPanelAdaptor(yawSetter1, viewPanel);
    yawSetter1.addChangeListener(MToP); 

    ModelToPanelAdaptor MToP2 = new ModelToPanelAdaptor(yawSetter2, viewPanel);
    yawSetter2.addChangeListener(MToP2); 
  }

//---------------------------------------

   public static void main(String args[]) {

   //-------------------------------------

      // the test bed
      test frame = new test("Yaw");


      frame.addWindowListener(new WindowAdapter() {
	public void windowClosing(WindowEvent e) {
	  System.exit(0);
	}
      });

      frame.setSize(800, 150);
      frame.setVisible(true);
   } // main()

} // class test

