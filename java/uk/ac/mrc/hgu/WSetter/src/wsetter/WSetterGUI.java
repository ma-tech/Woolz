package wsetter;

import wsetter.*;

import javax.swing.*;
import java.awt.*;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;

/**
 * Super class for WSetter.
 * <br>Separates out the GUI stuff for JBuilder compatibility
 * @author Nick Burton
 * @see WSetter
 */
public class WSetterGUI extends JPanel {

   // GUI stuff
   protected JPanel textPanel = new JPanel();
   protected JTextField textf = new JTextField();

   protected JPanel sliderPanel = new JPanel();
   protected JSlider _slider = new JSlider();

   protected JPanel labelSliderTextPanel = new JPanel();

   protected JPanel labelPanel = new JPanel();
   protected JLabel paramLabel = new JLabel("", SwingConstants.CENTER);

   /**
    * The nominal width of the WSetter.
    * This field is required as a work-around for a GUI problem.
    * The width of the WSetter is determined by Java Layout Manager, however
    * if it is set to say 800 and the window widened >> 800 it doesn't
    * draw itself beyond 800? unless you click on it.
    * If it is set to a number > screen resolution, say 2048, everything is OK
    */
   protected int _width = 2048;
   /**
    * The height of the WSetter in pixels.
    */
   protected int _height = 15;
   /**
    * The width of the WSetter text field in pixels.
    */
   protected int textWidth = 40;
   /**
    * The width of the WSetter label in pixels.
    */
   protected int labelWidth = 60;
   /**
    */
   protected int pad = 1;
   /**
    *
    */
   protected int _minWidth = 150; // not accessible
   protected int _minHeight = 10; // not accessible

   /**
    * The background colour of the WSetter text field and label.
    */
   protected Color bgc = new Color(230, 230, 230);
   /**
    * The background colour of the WSetter.
    */
   protected Color internalBgc = new Color(200, 200, 200); // not accessible

   protected String sliderLabel = "";
   // constants to define GUI shape
   protected int sliderH;
   protected int sliderW;
   //....................................

   /**
    *
    */
   public WSetterGUI() {
      initGUI();
   }

   /**
    * Initialises the WSetter GUI
    * @param	void
    * @return	void
    */
   protected void initGUI() {
      sliderH = _height;

      this.setLayout(new BorderLayout());
      this.setBackground(internalBgc);

      textf.setBackground(bgc);
      textf.setFont(new Font("default", Font.PLAIN, 11));
      textf.setPreferredSize(new Dimension(textWidth, sliderH+pad));
      textf.setMinimumSize(new Dimension(textWidth, sliderH+pad));
      textf.setMaximumSize(new Dimension(textWidth, sliderH+pad));

/*
      textPanel.setBackground(internalBgc);
      textPanel.setLayout(new BorderLayout(0, pad));
      textPanel.setPreferredSize(new Dimension(textWidth+pad, sliderH+pad));
      textPanel.add(textf, BorderLayout.NORTH);
*/

      sliderPanel.setBackground(internalBgc);
      sliderPanel.setLayout(new BorderLayout());
      sliderPanel.add(_slider, BorderLayout.CENTER);

      //paramLabel.setBackground(bgc);
      paramLabel.setFont(new Font("default", Font.PLAIN, 11));
      paramLabel.setPreferredSize(new Dimension(labelWidth, sliderH));
      paramLabel.setText(sliderLabel);

      labelPanel.setBackground(bgc);
      labelPanel.setLayout(new BorderLayout(pad, pad));
      labelPanel.setPreferredSize(new Dimension(labelWidth, sliderH));
      labelPanel.setMaximumSize(new Dimension(labelWidth, sliderH));
      labelPanel.setMinimumSize(new Dimension(labelWidth, sliderH));
      labelPanel.add(paramLabel, BorderLayout.NORTH);

      labelSliderTextPanel.setBackground(internalBgc);
      labelSliderTextPanel.setLayout(new BoxLayout(labelSliderTextPanel, BoxLayout.LINE_AXIS));
      labelSliderTextPanel.add(labelPanel);
      labelSliderTextPanel.add(sliderPanel);
      labelSliderTextPanel.add(textf);
      this.add(labelSliderTextPanel, BorderLayout.NORTH);

   } // initGUI

} // class WSetter
