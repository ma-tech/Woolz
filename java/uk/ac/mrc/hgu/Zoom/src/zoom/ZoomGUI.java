package zoom;

import zoom.*;

import javax.swing.*;
import java.awt.*;
import javax.swing.border.*;
import java.awt.event.*;
import javax.swing.event.*;
import java.io.*;
import java.net.*;

/**
 * GUI super class for Zoom.
 * <br>Separates out the GUI stuff for tidiness and JBuilder woossies
 * @author Nick Burton
 * @see Zoom
 */

public class ZoomGUI extends JComponent
                     implements Serializable {

   // GUI stuff

   protected JPanel topPanel = new JPanel();

   protected JPanel textPanel = new JPanel();
   protected JTextField textf = new JTextField();

   protected JPanel buttonPanel = new JPanel();
   protected JPanel buttonPanelIn = new JPanel();
   protected JPanel buttonPanelOut = new JPanel();


   protected ImageIcon _inIcon = null;
   protected ImageIcon _outIcon = null;

   protected JButton _inButton = null;
   protected JButton _outButton = null;

   protected JPanel labelTextPanel = new JPanel();
   protected JPanel labelTextButtonPanel = new JPanel();

   protected JPanel labelPanel = new JPanel();
   protected JLabel paramLabel = new JLabel("", SwingConstants.CENTER);

   // constants to define GUI shape
   //protected int _height = 15;
   //protected int _width = 150;
   protected int _buttonH = 15;
   protected int _buttonW = 25;
   protected int _textWidth = 50;
   protected int _labelWidth = 50;
   protected int _pad = 1;
   protected int _minWidth = 150; // not accessible
   protected int _maxWidth = 300; // not accessible
   protected int _minHeight = 15; // not accessible
   protected int _maxHeight = 30; // not accessible

   protected Color _bgc = new Color(230, 230, 230);
   protected Color _internalBgc = new Color(200, 200, 200); // not accessible

   protected String _zoomLabel = "zoom %";
   //....................................

   /**
    * Constructor
    */
   public ZoomGUI() {
      initGUI();
   }

   /**
    * Initialises the GUI
    * @param	void
    * @return	void
    */
   protected void initGUI() {

      this.setLayout(new BorderLayout());
      this.setBackground(_internalBgc);

      textf.setBackground(_bgc);
      textf.setPreferredSize(new Dimension(_textWidth, _buttonH+_pad));

      textPanel.setBackground(_internalBgc);
      textPanel.setLayout(new BorderLayout(0, _pad));
      textPanel.setPreferredSize(new Dimension(_textWidth+_pad, _buttonH+_pad));
      textPanel.add(textf, BorderLayout.NORTH);

      _outIcon = null;
      _inIcon = null;
      //URL urlIn = ZoomGUI.class.getResource("images/in.gif");
      //URL urlOut = ZoomGUI.class.getResource("images/out.gif");

      // Get current classloader
      // this is needed to retrieve resources from a jar file (ie from web start deployed application)
      ClassLoader cl = ZoomGUI.class.getClassLoader();
      String imgPath = "images/";
      _outIcon = new ImageIcon(cl.getResource(imgPath + "out.gif"));
      _inIcon = new ImageIcon(cl.getResource(imgPath + "in.gif"));


      //if(urlIn != null) _inIcon = new ImageIcon(urlIn);
      //if(urlOut != null) _outIcon = new ImageIcon(urlOut);
      if(_inIcon != null) {
	 _inButton = new JButton(_inIcon);
      } else {
	 _inButton = new JButton("+");
      }
      if(_outIcon != null) {
	 _outButton = new JButton(_outIcon);
      } else {
	 _outButton = new JButton("-");
      }

      buttonPanelIn.setBackground(_internalBgc);
      buttonPanelIn.setLayout(new BorderLayout());
      buttonPanelIn.setPreferredSize(
                      new Dimension(_buttonW+2*_pad, _buttonH+_pad));
      buttonPanelIn.add(_inButton, BorderLayout.CENTER);

      buttonPanelOut.setBackground(_internalBgc);
      buttonPanelOut.setLayout(new BorderLayout());
      buttonPanelOut.setPreferredSize(
                      new Dimension(_buttonW+2*_pad, _buttonH+_pad));
      buttonPanelOut.add(_outButton, BorderLayout.CENTER);

      buttonPanel.setBackground(_internalBgc);
      buttonPanel.setLayout(new BorderLayout());
      buttonPanel.setPreferredSize(new Dimension(
                      2*_buttonW+6*_pad, _buttonH+_pad));
      buttonPanel.add(buttonPanelOut, BorderLayout.WEST);
      buttonPanel.add(buttonPanelIn, BorderLayout.EAST);

      paramLabel.setFont(new Font("default", Font.PLAIN, 11));
      paramLabel.setBackground(_internalBgc);
      paramLabel.setPreferredSize(new Dimension(_labelWidth, _buttonH));
      paramLabel.setText(_zoomLabel);

      labelPanel.setBackground(_bgc);
      labelPanel.setLayout(new BorderLayout(_pad, _pad));
      labelPanel.setPreferredSize(new Dimension(_labelWidth, _buttonH));
      labelPanel.add(paramLabel, BorderLayout.NORTH);

      labelTextPanel.setBackground(_internalBgc);
      labelTextPanel.setLayout(new BorderLayout(_pad,_pad));
      labelTextPanel.add(labelPanel, BorderLayout.WEST);
      labelTextPanel.add(textPanel, BorderLayout.EAST);

      labelTextButtonPanel.setBackground(_internalBgc);
      labelTextButtonPanel.setLayout(new BorderLayout());
      labelTextButtonPanel.setPreferredSize(
         new Dimension(_labelWidth+_textWidth+(2*_buttonW)+(6*_pad), _buttonH+_pad));
      labelTextButtonPanel.add(labelTextPanel, BorderLayout.WEST);
      labelTextButtonPanel.add(buttonPanel, BorderLayout.EAST);

      topPanel.setBackground(_bgc);
      topPanel.setLayout(new BorderLayout(_pad, _pad));
/*
      topPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createEmptyBorder(1,1,1,1),
                BorderFactory.createEtchedBorder(EtchedBorder.RAISED)));
*/
      topPanel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.RAISED));
      topPanel.add(labelTextButtonPanel, BorderLayout.WEST);

      this.add(topPanel, BorderLayout.WEST);

   } // initGUI

//.....................................
// the following are required for setting the bean properties
//.....................................
   public void setLabelWidth(int w) {
      _labelWidth = w;
      if(labelPanel != null)
         labelPanel.setPreferredSize(new Dimension(_labelWidth, _buttonH));
      if(labelTextButtonPanel != null)
         labelTextButtonPanel.setPreferredSize(
                  new Dimension(_labelWidth+_textWidth+(2*_buttonW)+(6*_pad), _buttonH+_pad));
   }
   public void setTextWidth(int w) {
      _textWidth = w;
      if(textf != null)
         textf.setPreferredSize(new Dimension(_textWidth, _buttonH+_pad));
      if(textPanel != null)
         textPanel.setPreferredSize(new Dimension(_textWidth+_pad, _buttonH+_pad));
      if(labelTextButtonPanel != null)
         labelTextButtonPanel.setPreferredSize(
                  new Dimension(_labelWidth+_textWidth+(2*_buttonW)+(6*_pad), _buttonH+4*_pad));
   }
   public void setBgc(Color col) {
      _bgc = col;
      if(labelPanel != null) labelPanel.setBackground(_bgc);
      if(topPanel != null) topPanel.setBackground(_bgc);
      if(textf != null) textf.setBackground(_bgc);
   }
   public void setZoomLabel(String labl) {
      _zoomLabel = labl;
      if(paramLabel != null) paramLabel.setText(_zoomLabel);
   }
//.....................................

} // class Zoom
