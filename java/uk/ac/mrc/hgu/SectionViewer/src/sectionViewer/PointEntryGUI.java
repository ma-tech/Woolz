package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;

/**
 *   GUI for the PointEntry class
 */
public class PointEntryGUI extends JFrame{

  protected JPanel topPanel;
  protected JPanel msgPanel;
  protected JPanel buttonPanel;
  protected JPanel XPanel;
  protected JPanel YPanel;
  protected JPanel ZPanel;
  protected JPanel XPanel1;
  protected JPanel YPanel1;
  protected JPanel ZPanel1;

  protected Box XYZBox;
  protected Box controlBox;

  protected JPanel tfPanel1;
  protected JPanel tfPanel2;
  protected JPanel tfPanel3;

  protected JPanel lblPanel1;
  protected JPanel lblPanel2;
  protected JPanel lblPanel3;

  protected JLabel lblMsg;
  protected JLabel lblX;
  protected JLabel lblY;
  protected JLabel lblZ;

  protected JTextField tfX;
  protected JTextField tfY;
  protected JTextField tfZ;

  protected JButton _okButton;
  protected JButton _cancelButton;
  protected String _msg = "";
  protected int _gap;
//-------------------------------------------------------------
  /**
   *   Creates a PointEntry GUI with the given title.
   *   @param msg the title of the PointEntry GUI.
   */
  public PointEntryGUI(String msg) {
    super("PointEntry");

    topPanel = (JPanel)this.getContentPane();
    topPanel.setPreferredSize(new Dimension(200,150));
    _gap = 3;
    _msg = msg;
    makeGUI();
  }

//-------------------------------------------------------------
  /**
   *   Most of the work of building the GUI is done here.
   */
  public void makeGUI() {
    msgPanel = new JPanel();
    msgPanel.setLayout(new FlowLayout());
    msgPanel.setBackground(new Color(220,220,220));
    msgPanel.setPreferredSize(new Dimension(0,30));
    msgPanel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED));

    controlBox = new Box(BoxLayout.Y_AXIS);

    buttonPanel = new JPanel();
    buttonPanel.setLayout(new GridLayout(1,0,_gap,_gap));
    buttonPanel.setPreferredSize(new Dimension(0,30));

    //-----------------------------
    lblPanel1 = new JPanel();
    lblPanel1.setLayout(new BorderLayout(_gap*3, _gap));

    tfPanel1 = new JPanel();
    tfPanel1.setLayout(new BorderLayout(_gap*3, _gap));

    XPanel = new JPanel();
    XPanel.setLayout(new FlowLayout());
    XPanel.setBackground(new Color(200,200,200));
    //-----------------------------
    lblPanel2 = new JPanel();
    lblPanel2.setLayout(new BorderLayout(_gap*3, _gap));

    tfPanel2 = new JPanel();
    tfPanel2.setLayout(new BorderLayout(_gap*3, _gap));

    YPanel = new JPanel();
    YPanel.setLayout(new FlowLayout());
    YPanel.setBackground(new Color(200,200,200));
    //-----------------------------
    lblPanel3 = new JPanel();
    lblPanel3.setLayout(new BorderLayout(_gap*3, _gap));

    tfPanel3 = new JPanel();
    tfPanel3.setLayout(new BorderLayout(_gap*3, _gap));

    ZPanel = new JPanel();
    ZPanel.setLayout(new FlowLayout());
    ZPanel.setBackground(new Color(200,200,200));
    //-----------------------------

    XYZBox = new Box(BoxLayout.Y_AXIS);
    //-----------------------------

    lblMsg = new JLabel(_msg);
    lblX = new JLabel("X");
    lblY = new JLabel("Y");
    lblZ = new JLabel("Z");
    //-----------------------------

    tfX = new JTextField("0",10);
    tfY = new JTextField("0",10);
    tfZ = new JTextField("0",10);
    //-----------------------------

    _okButton = new JButton("OK");
    _cancelButton = new JButton("Cancel");
    _okButton.setPreferredSize(new Dimension(0,20));
    _cancelButton.setPreferredSize(new Dimension(0,20));
    //-----------------------------
    msgPanel.add(lblMsg);
    //.............................
    //.............................
    buttonPanel.add(_okButton);
    buttonPanel.add(_cancelButton);
    //.............................
    lblPanel1.add(lblX, BorderLayout.WEST);
    tfPanel1.add(tfX, BorderLayout.EAST);
    XPanel.add(lblPanel1);
    XPanel.add(tfPanel1);
    //.............................
    lblPanel2.add(lblY, BorderLayout.WEST);
    tfPanel2.add(tfY, BorderLayout.EAST);
    YPanel.add(lblPanel2);
    YPanel.add(tfPanel2);
    //.............................
    lblPanel3.add(lblZ, BorderLayout.WEST);
    tfPanel3.add(tfZ, BorderLayout.EAST);
    ZPanel.add(lblPanel3);
    ZPanel.add(tfPanel3);
    //.............................
    XYZBox.add(XPanel);
    XYZBox.add(YPanel);
    XYZBox.add(ZPanel);
    //.............................
    controlBox.add(msgPanel);
    controlBox.add(XYZBox);
    controlBox.add(buttonPanel);
    //.............................
    topPanel.add(controlBox);
  }
//-------------------------------------------------------------
  /**
  *   Abstract Event Handler (to be implemented in PointEntry class).
  */
  abstract class ButtonHandler implements ActionListener {}

  /**
  *   Abstract Event Handler (to be implemented in PointEntry class).
  */
  abstract class TextFieldHandler implements ActionListener {}
//-------------------------------------------------------------

} // class PointEntryGUI
