package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

import wsetter.*;
import zoom.*;

import hguUntil.*;

public class SectionViewerGUI extends JPanel {

  private final boolean _debug = false;

  // these need to be visible outside of guiInit()
  int totalH = 500;
  int totalW = 400;
  int minW = 350;
  int minH = 300;
  int bord = 2;
  int pad = 1;
  int hgap = 2;
  int vgap = 2;
  int titleH;
  int pyrH;
  int rotH;
  int transH;

  Font titleFont;
  Font menuFont;
  Font feedbackFont;
  Font anatomyFont;

  JPanel topPanel = new JPanel();

  JPanel titlePanel = new JPanel();
  JPanel titlePanel_1 = new JPanel();
  JPanel titlePanel_2 = new JPanel();
  JPanel titlePanel_3 = new JPanel();

  ScrollableTextField titleText = new ScrollableTextField();
  JPanel menuPanel = new JPanel();
  JPanel titleMenuPanel = new JPanel();

  JPanel feedbackPanel = new JPanel();
  JPanel imagePanel = new JPanel();
  JPanel _bigPanel = new JPanel();
  JScrollPane _imageScrollPane = new JScrollPane();
  JPanel feedbackImagePanel = new JPanel();

  Zoom zoomSetter = new Zoom();
  WSetter distSetter = new WSetter();
  WSetter pitchSetter = new WSetter();
  WSetter yawSetter = new WSetter();
  WSetter rollSetter = new WSetter();
  WSetter rotSetter = new WSetter();

  JPanel zoomPanel = new JPanel();

  JPanel zoomControlPanel = new JPanel();
  JPanel spacerPanel_R = new JPanel();
  JPanel basicControlPanel = new JPanel();
  JPanel pitchYawRollPanel = new JPanel();
  JPanel pitchyawControlPanel = new JPanel();
  JPanel rollControlPanel = new JPanel();

  JPanel fixedPointControlPanel = new JPanel();
  JPanel userDefinedRotPanel = new JPanel();
  JPanel fixedPointUserRotPanel = new JPanel();

  JPanel menuPanel_2 = new JPanel();
  JMenuBar _menubar = new JMenuBar();
  //...................

  JPanel permanentPanel = new JPanel();
  JPanel transientPanel = new JPanel();

  //...................
  String fileMenuStr = "file";
  String fileMenu_1str = "save image";
  String fileMenu_2str = "save view settings";
  String fileMenu_3str = "load view settings";
  String fileMenu_4str = "close view";
  String fileMenu_5str = "save anatomy image";

  JMenu fileMenu = new JMenu(fileMenuStr);
  JMenuItem fileMenu_1 = new JMenuItem(fileMenu_1str);
  JMenuItem fileMenu_2 = new JMenuItem(fileMenu_2str);
  JMenuItem fileMenu_3 = new JMenuItem(fileMenu_3str);
  JMenuItem fileMenu_4 = new JMenuItem(fileMenu_4str);
  JMenuItem fileMenu_5 = new JMenuItem(fileMenu_5str);
  //...................

  String controlMenuStr = "control";
  String controlMenu_1str = "show rotation controls";
  String controlMenu_2str = "show user-defined rotation control";
  String controlMenu_3str = "reset controls";

  String viewMenuStr = "view";
  String viewMenu_1str = "show cursor feedback";
  String viewMenu_2str = "show intersection of views";
  String viewMenu_3str = "enable 'mouse-click' anatomy";

  JMenu controlMenu = new JMenu(controlMenuStr);
  boolean controlMenu_1state = false;
  JCheckBoxMenuItem controlMenu_1 = new
               JCheckBoxMenuItem(controlMenu_1str,
	                         controlMenu_1state);
  boolean controlMenu_2state = false;
  JCheckBoxMenuItem controlMenu_2 = new
               JCheckBoxMenuItem(controlMenu_2str,
	                         controlMenu_2state);

  JMenuItem controlMenu_3 = new JMenuItem(controlMenu_3str);

  JMenu viewMenu = new JMenu(viewMenuStr);
  boolean viewMenu_1state = true;
  JCheckBoxMenuItem viewMenu_1 = new
               JCheckBoxMenuItem(viewMenu_1str,
                                 viewMenu_1state);

  boolean viewMenu_2state = false;
  JCheckBoxMenuItem viewMenu_2 = new
               JCheckBoxMenuItem(viewMenu_2str,
                                 viewMenu_2state);

  boolean viewMenu_3state = false;
  JCheckBoxMenuItem viewMenu_3 = new
               JCheckBoxMenuItem(viewMenu_3str,
                                 viewMenu_3state);


  //...................
  String viewModeMenuStr = "view_mode";
  String viewModeMenu_1str = "up_is_up";
  String viewModeMenu_2str = "absolute";
  //String viewModeMenu_3str = "help";

  JMenu viewModeMenu = new JMenu(viewModeMenuStr);
  ButtonGroup viewModeGroup = new ButtonGroup();
  JRadioButtonMenuItem viewModeMenu_1 = new
               JRadioButtonMenuItem(viewModeMenu_1str);
  JRadioButtonMenuItem viewModeMenu_2 = new
               JRadioButtonMenuItem(viewModeMenu_2str);
  //JMenuItem viewModeMenu_3 = new JMenuItem(viewModeMenu_3str);
  //...................
  String thresholdMenuStr = "threshold";
  String thresholdMenu_1str = "enable constraint definition";
  String thresholdMenu_2str = "enable thresholding";
  String thresholdMenu_3str = "show threshold constraint";
  String thresholdMenu_4str = "remove threshold constraint";

  boolean thresholdMenu_1state = false;
  boolean thresholdMenu_2state = false;
  boolean thresholdMenu_3state = false;
  JMenu thresholdMenu = new JMenu(thresholdMenuStr);
  JCheckBoxMenuItem thresholdMenu_1 = new
               JCheckBoxMenuItem(thresholdMenu_1str, thresholdMenu_1state);
  JCheckBoxMenuItem thresholdMenu_2 = new
               JCheckBoxMenuItem(thresholdMenu_2str, thresholdMenu_2state);
  JCheckBoxMenuItem thresholdMenu_3 = new
               JCheckBoxMenuItem(thresholdMenu_3str, thresholdMenu_3state);
  JMenuItem thresholdMenu_4 = new JMenuItem(thresholdMenu_4str);
  //...................
  String fixedPointMenuStr = "fixed_point";
  String fixedPointMenu_1str = "change fixed point";
  String fixedPointMenu_1_1str = "using left mouse button";
  String fixedPointMenu_1_2str = "typing coordinates";
  String fixedPointMenu_2str = "reset fixed point";
  String fixedPointMenu_3str = "show fixed point";
  String fixedPointMenu_4str = "define rotation axis";
  String fixedPointMenu_4_1str = "using left mouse button";
  String fixedPointMenu_4_2str = "typing coordinates";
  String fixedPointMenu_5str = "display rotation axis";
  //String fixedPointMenu_6str = "help";

  JMenu fixedPointMenu = new JMenu(fixedPointMenuStr);
  boolean fixedPointMenu_3state = false;
  JMenu fixedPointMenu_1 = new JMenu(fixedPointMenu_1str);
  JMenuItem fixedPointMenu_1_1 = new JMenuItem(fixedPointMenu_1_1str);
  JMenuItem fixedPointMenu_1_2 = new JMenuItem(fixedPointMenu_1_2str);
  JMenuItem fixedPointMenu_2 = new JMenuItem(fixedPointMenu_2str);
  JCheckBoxMenuItem fixedPointMenu_3 = new
               JCheckBoxMenuItem(fixedPointMenu_3str,
	                         fixedPointMenu_3state);
  JMenu fixedPointMenu_4 = new JMenu(fixedPointMenu_4str);
  JMenuItem fixedPointMenu_4_1 = new JMenuItem(fixedPointMenu_4_1str);
  JMenuItem fixedPointMenu_4_2 = new JMenuItem(fixedPointMenu_4_2str);
  JCheckBoxMenuItem fixedPointMenu_5 = new
               JCheckBoxMenuItem(fixedPointMenu_5str, true);
  //JMenuItem fixedPointMenu_6 = new JMenuItem(fixedPointMenu_6str);
  //...................

  String helpMenuStr = "help";
  String helpMenu_1str = "contents";
  String helpMenu_2str = "index";
  String helpMenu_3str = "search";
  JMenuItem helpMenu_1 = new JMenuItem(helpMenu_1str);
  JMenuItem helpMenu_2 = new JMenuItem(helpMenu_2str);
  JMenuItem helpMenu_3 = new JMenuItem(helpMenu_3str);

  JMenu helpMenu = new JMenu(helpMenuStr);
  //...................
  JPanel feedbackPanel_1 = new JPanel();
  JPanel feedbackPanel_1a = new JPanel();
  JPanel feedbackPanel_2 = new JPanel();
  JPanel feedbackPanel_3 = new JPanel();
  JPanel feedbackPanel_4 = new JPanel();
  JPanel feedbackPanel_5 = new JPanel();
  JPanel feedbackPanel_6 = new JPanel();
  JPanel feedbackPanel_7 = new JPanel();
  JPanel fbAnatPanel = new JPanel();

  JTextField xyzTextField = new JTextField();
  JTextField valueTextField = new JTextField();
  ScrollableTextField anatomyTextField = new ScrollableTextField();

  String titleString = "Ref:";

  Color titleCol = new Color(200, 220, 220);
  Color fbCol = new Color(220, 220, 200);

  JButton secColorClt = new JButton("");

  static int nSV = 0;

  private Color[] planeColor = new Color[] {Color.red, Color.yellow,
      Color.blue, Color.pink, Color.green, Color.cyan};
  //=========================================================
  // constructor
  //=========================================================
  public SectionViewerGUI() {

    if(_debug == true) System.out.println("enter SectionViewerGUI");

    try {
      guiInit();
    }
    catch(Exception e) {
      System.out.println("SectionViewerGUI");
      e.printStackTrace();
    }

    if(_debug) System.out.println("exit SectionViewer");
  }

//========================================
  private void guiInit() throws Exception {

    if(_debug) System.out.println("enter guiInit");

    // constants to define GUI shape
    int menuH;

    int titleW = 200;
    int fbH;
    int xyzW = 105;
    int xyzH;
    int valW = 65;
    int valH;
    int anlabW = 60;
    int antextW = 150;
    int anatW = anlabW+antextW+hgap;
    int anatH;

    int imgH = 200;
    int fbimgH;
    int distH = 25;

    int permH;
    //......................................
    int pyH = 40;
    int rH = 20;
    pyrH = pyH+rH;
    rotH = 25;
    transH = pyrH+vgap+rotH;

//----------------------------------------------------------
    fileMenu.add(fileMenu_1);
    fileMenu.add(fileMenu_5);
    fileMenu.addSeparator();
    fileMenu.add(fileMenu_2);
    fileMenu.add(fileMenu_3);
    fileMenu.addSeparator();
    fileMenu.add(fileMenu_4);
    fileMenu_1.setEnabled(true);
    fileMenu_2.setEnabled(true);
    fileMenu_3.setEnabled(true);
    fileMenu_4.setEnabled(true);
 //...............................
    controlMenu.add(controlMenu_1);
    controlMenu.add(controlMenu_2);
    controlMenu.add(controlMenu_3);
    controlMenu.addSeparator();

    viewMenu.add(viewMenu_1);
    viewMenu.add(viewMenu_2);
    viewMenu.add(viewMenu_3);
    viewMenu_3.setEnabled(false);

    controlMenu.add(viewModeMenu_1);
    controlMenu.add(viewModeMenu_2);
    controlMenu.addSeparator();

 //...............................
    viewModeGroup.add(viewModeMenu_1);
    viewModeGroup.add(viewModeMenu_2);
    viewModeMenu_1.setSelected(true);

    controlMenu.add(fixedPointMenu_1);
    controlMenu.add(fixedPointMenu_2);
    controlMenu.add(fixedPointMenu_3);
    controlMenu.addSeparator();
    controlMenu.add(fixedPointMenu_4);
    controlMenu.add(fixedPointMenu_5);

 //...............................
    fixedPointMenu_1.add(fixedPointMenu_1_1);
    fixedPointMenu_1.add(fixedPointMenu_1_2);
    fixedPointMenu_4.add(fixedPointMenu_4_1);
    fixedPointMenu_4.add(fixedPointMenu_4_2);
    fixedPointMenu.setEnabled(true);
    fixedPointMenu_1_1.setEnabled(true);
    fixedPointMenu_1_2.setEnabled(true);
    fixedPointMenu_2.setEnabled(true);
    fixedPointMenu_3.setEnabled(true);
    fixedPointMenu_4.setEnabled(false);
    fixedPointMenu_4_1.setEnabled(false);
    fixedPointMenu_4_2.setEnabled(false);
    fixedPointMenu_5.setEnabled(false);
 //...............................
    thresholdMenu.add(thresholdMenu_1);
    thresholdMenu.add(thresholdMenu_2);
    thresholdMenu.add(thresholdMenu_3);
    thresholdMenu.add(thresholdMenu_4);
    thresholdMenu.setEnabled(true);
    thresholdMenu_1.setEnabled(true);
    thresholdMenu_2.setEnabled(true);
    thresholdMenu_3.setEnabled(true);
    thresholdMenu_4.setEnabled(true);
 //...............................
    helpMenu.add(helpMenu_1);
    helpMenu.add(helpMenu_2);
    helpMenu.add(helpMenu_3);
    helpMenu_1.setEnabled(true);
    helpMenu_2.setEnabled(true);
    helpMenu_3.setEnabled(true);
    //----------------------------------------------------------
    _menubar.add(fileMenu);
    _menubar.add(controlMenu);
    _menubar.add(viewMenu);
    //_menubar.add(thresholdMenu);
    _menubar.add(helpMenu);
    //----------------------------------------------------------
    titleFont = new Font("Helvetica", Font.PLAIN, 12);
    titleH = titleText.getFontMetrics(titleFont).getHeight() + bord;
    //......................................
    menuFont = new Font("Helvetica", Font.PLAIN, 12);
    menuH = _menubar.getFontMetrics(menuFont).getHeight() + bord;
    //......................................
    feedbackFont = new Font("Helvetica", Font.PLAIN, 12);
    anatomyFont = new Font("Helvetica", Font.ITALIC, 12);
    fbH = titleText.getFontMetrics(feedbackFont).getHeight() + bord;
    //......................................
    xyzH = fbH+2*bord;
    valH = fbH+2*bord;
    anatH = fbH+2*bord;
    fbimgH = fbH+imgH+vgap;
    permH = menuH+vgap+fbimgH+vgap+distH;
    //----------------------------------------------------------

    titleText.setEditable(false);
    titleText.setHorizontalAlignment(JTextField.CENTER);
    titleText.setPreferredSize(new Dimension(30,titleH));
    titleText.setFont(titleFont);
    titleText.setText("no grey level file");
    titleText.setBackground(titleCol);

    xyzTextField.setFont(feedbackFont);
    xyzTextField.setEditable(false);
    xyzTextField.setBackground(fbCol);

    valueTextField.setFont(feedbackFont);
    valueTextField.setEditable(false);
    valueTextField.setBackground(fbCol);

    anatomyTextField.setFont(feedbackFont);
    anatomyTextField.setEditable(false);
    anatomyTextField.setBackground(fbCol);

    titlePanel_1.setPreferredSize(new Dimension(50,titleH+4*bord));
    titlePanel_1.setLayout(new BorderLayout(hgap, vgap));
    titlePanel_1.setBorder(BorderFactory.createEmptyBorder(1,1,1,1));
    titlePanel_1.add(titleText, BorderLayout.CENTER);

    titlePanel_2.setPreferredSize(new Dimension(50, 0));
    titlePanel_3.setPreferredSize(new Dimension(50, 0));

    titlePanel.setPreferredSize(new Dimension(50,titleH+2*bord));
    titlePanel.setLayout(new BorderLayout(hgap, vgap));
    titlePanel.add(titlePanel_1, BorderLayout.CENTER);
    titlePanel.add(titlePanel_2, BorderLayout.WEST);
    titlePanel.add(titlePanel_3, BorderLayout.EAST);
    //...............................

    _menubar.setPreferredSize(new Dimension(0,menuH+2*bord));
    menuPanel.setPreferredSize(new Dimension(totalW, menuH+4*bord));
    menuPanel.setLayout(new BorderLayout(hgap, vgap));
    menuPanel_2.setPreferredSize(new Dimension(totalW, menuH+4*bord));
    menuPanel_2.setLayout(new BorderLayout(hgap, vgap));

    //...............................
    menuPanel_2.add(_menubar, BorderLayout.CENTER);

    menuPanel.add(menuPanel_2, BorderLayout.CENTER);

    //...............................
    titleMenuPanel.setLayout(new BorderLayout());
    titleMenuPanel.add(menuPanel, BorderLayout.SOUTH);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // the top level feedback panel
    feedbackPanel.setPreferredSize(new Dimension(totalW, fbH+4*bord));
    feedbackPanel.setLayout(new BorderLayout(hgap, vgap));
    feedbackPanel.setBorder(BorderFactory.createEtchedBorder(EtchedBorder.LOWERED));
    //...............................

    // panel for x,y info which fixes its width
    feedbackPanel_1.setPreferredSize(new Dimension(100, fbH + 2 * bord));
    feedbackPanel_1.setLayout(new BorderLayout(hgap, vgap));
    feedbackPanel_1.add(xyzTextField, BorderLayout.CENTER);

    // buffer between xyz and grey value feedback
    feedbackPanel_1a.setPreferredSize(new Dimension(5, fbH+2*bord));
    // panel for grey value info which fixes its width
    feedbackPanel_2.setPreferredSize(new Dimension(35, fbH+2*bord));
    feedbackPanel_2.setLayout(new BorderLayout(hgap, vgap));
    feedbackPanel_2.add(feedbackPanel_1a, BorderLayout.WEST);
    feedbackPanel_2.add(valueTextField, BorderLayout.CENTER);

    // combine position & value info
    feedbackPanel_3.setPreferredSize(new Dimension(130, fbH+2*bord));
    feedbackPanel_3.setLayout(new BorderLayout(hgap, vgap));
    feedbackPanel_3.add(feedbackPanel_1, BorderLayout.WEST);
    feedbackPanel_3.add(feedbackPanel_2, BorderLayout.EAST);

    //...............................

    // panel for anatomy text which lets it expand
    fbAnatPanel.setPreferredSize(new Dimension(0, fbH+2*bord));
    fbAnatPanel.setLayout(new BorderLayout(hgap, vgap));
    fbAnatPanel.setBorder(BorderFactory.createEmptyBorder(1,1,1,1));
    fbAnatPanel.add(anatomyTextField, BorderLayout.CENTER);
    //...............................

    // set up the colour chooser button
    secColorClt.setPreferredSize(new Dimension(titleH, titleH));
    secColorClt.setBorder(null);
    secColorClt.setFocusPainted(false);
    try {
      secColorClt.setBackground(new Color(
             planeColor[nSV % planeColor.length].getRed(),
             planeColor[nSV % planeColor.length].getGreen(),
             planeColor[nSV % planeColor.length].getBlue()));
      nSV++;
    }
    catch (Exception exp) {
      nSV = 0;
    }

    feedbackPanel.add(feedbackPanel_3, BorderLayout.WEST);
    feedbackPanel.add(fbAnatPanel, BorderLayout.CENTER);
    feedbackPanel.add(secColorClt, BorderLayout.EAST);
    //------------------------------
    // if this panel is bigger than the _imageScrollPane,
    // scroll bars should appear
    //...............................
    _bigPanel.setLayout(new BorderLayout());
    //...............................

    _imageScrollPane.getViewport().add(_bigPanel, BorderLayout.CENTER);
    _imageScrollPane.getViewport().setScrollMode(JViewport.BLIT_SCROLL_MODE);

    imagePanel.setLayout(new BorderLayout(hgap, vgap));
    imagePanel.setBorder(BorderFactory.createCompoundBorder(
                 BorderFactory.createEmptyBorder(1,1,1,1),
		 BorderFactory.createEtchedBorder(EtchedBorder.LOWERED)));
    imagePanel.add(_imageScrollPane, BorderLayout.CENTER);

    //...............................

    feedbackImagePanel.setLayout(new BorderLayout(hgap, vgap));
    feedbackImagePanel.add(imagePanel, BorderLayout.CENTER);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    zoomSetter.setBgc(new Color(230, 230, 230));
    zoomSetter.setMin(20);
    zoomSetter.setMax(1000);
    zoomSetter.setValue(100);
    zoomSetter.setInc(50);

    distSetter.setBgc(new Color(255, 255, 230));
    distSetter.setLabelWidth(60);
    distSetter.setTextWidth(40);
    distSetter.setSliderLabel("dist");

    zoomControlPanel.setLayout(new BorderLayout(hgap, vgap));
    zoomControlPanel.add(zoomSetter, BorderLayout.WEST);
    zoomControlPanel.add(spacerPanel_R, BorderLayout.CENTER);

    basicControlPanel.setPreferredSize(new Dimension(totalW, 2*distH));
    basicControlPanel.setLayout(new BorderLayout(hgap, vgap));
    basicControlPanel.setBorder(BorderFactory.createCompoundBorder(
                 BorderFactory.createEmptyBorder(3,1,3,1),
		             BorderFactory.createEtchedBorder(EtchedBorder.LOWERED)));
    basicControlPanel.add(zoomSetter, BorderLayout.NORTH);
    basicControlPanel.add(distSetter, BorderLayout.SOUTH);

    pitchSetter.setBgc(new Color(255, 230, 230));
    pitchSetter.setLabelWidth(60);
    pitchSetter.setTextWidth(40);
    pitchSetter.setSliderLabel("pitch");
    pitchSetter.setMin(0.0);
    pitchSetter.setMax(180.0);
    pitchSetter.setValue(0.0);

    yawSetter.setBgc(new Color(230, 255, 230));
    yawSetter.setLabelWidth(60);
    yawSetter.setTextWidth(40);
    yawSetter.setSliderLabel("yaw");
    yawSetter.setMin(0.0);
    yawSetter.setMax(180.0);
    yawSetter.setValue(0.0);

    pitchyawControlPanel.setPreferredSize(new Dimension(totalW, pyH));
    pitchyawControlPanel.setMinimumSize(new Dimension(totalW, pyH));
    pitchyawControlPanel.setLayout(new BorderLayout(hgap, vgap));
    pitchyawControlPanel.add(pitchSetter, BorderLayout.NORTH);
    pitchyawControlPanel.add(yawSetter, BorderLayout.SOUTH);

    rollSetter.setBgc(new Color(230, 230, 255));
    rollSetter.setLabelWidth(60);
    rollSetter.setTextWidth(40);
    rollSetter.setSliderLabel("roll");
    rollSetter.setMin(0.0);
    rollSetter.setMax(180.0);
    rollSetter.setValue(0.0);
    rollSetter.setSliderEnabled(false);

    rollControlPanel.setPreferredSize(new Dimension(totalW, rH));
    rollControlPanel.setMinimumSize(new Dimension(totalW, rH));
    rollControlPanel.setLayout(new BorderLayout(hgap, vgap));
    rollControlPanel.add(rollSetter, BorderLayout.SOUTH);

    pitchYawRollPanel.setPreferredSize(new Dimension(totalW, pyrH));
    pitchYawRollPanel.setMinimumSize(new Dimension(totalW, pyrH));
    pitchYawRollPanel.setLayout(new BorderLayout(hgap, vgap));
    pitchYawRollPanel.add(pitchyawControlPanel, BorderLayout.NORTH);
    pitchYawRollPanel.add(rollControlPanel, BorderLayout.CENTER);

    rotSetter.setBgc(new Color(230, 255, 255));
    rotSetter.setLabelWidth(60);
    rotSetter.setTextWidth(40);
    rotSetter.setSliderLabel("fixed rotn");
    rotSetter.setMin(0.0);
    rotSetter.setMax(360.0);
    rotSetter.setValue(0.0);
    rotSetter.setSliderEnabled(false);

    userDefinedRotPanel.setPreferredSize(new Dimension(totalW, rotH));
    userDefinedRotPanel.setMinimumSize(new Dimension(totalW, rotH));
    userDefinedRotPanel.setLayout(new BorderLayout(hgap, vgap));
    userDefinedRotPanel.setBorder(BorderFactory.createCompoundBorder(
                 BorderFactory.createEmptyBorder(3,1,3,1),
		             BorderFactory.createEtchedBorder(EtchedBorder.LOWERED)));
    userDefinedRotPanel.add(rotSetter, BorderLayout.CENTER);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    permanentPanel.setLayout(new BorderLayout());
    permanentPanel.add(titleMenuPanel, BorderLayout.NORTH);
    permanentPanel.add(feedbackImagePanel, BorderLayout.CENTER);
    permanentPanel.add(basicControlPanel, BorderLayout.SOUTH);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    transientPanel.setLayout(new BorderLayout());
    transientPanel.setPreferredSize(new Dimension(totalW, pad));
    transientPanel.setMinimumSize(new Dimension(totalW, pad));
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    topPanel.setLayout(new BorderLayout());
    topPanel.setPreferredSize(new Dimension(totalW+pad, totalH+pad));
    topPanel.setMinimumSize(new Dimension(totalW+pad, totalH+pad));
    topPanel.add(permanentPanel, BorderLayout.CENTER);
    topPanel.add(transientPanel, BorderLayout.SOUTH);

    //----------------------------------------------------------
    this.setLayout(new BorderLayout());
    this.add(topPanel, BorderLayout.CENTER);
    this.setPreferredSize(new Dimension(totalW+pad, totalH+pad));
    this.setBorder(BorderFactory.createCompoundBorder(
             BorderFactory.createEmptyBorder(3,1,3,1),
             BorderFactory.createEtchedBorder(EtchedBorder.LOWERED)));
    this.setVisible(true);
    //----------------------------------------------------------

    if(_debug) System.out.println("exit guiInit");
  } // guiInit()


//----------------------------------------------------------

//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
  public abstract class fileMenuHandler implements ActionListener {
  }
//---------------------------------------
  public abstract class controlMenuHandler implements ActionListener {
  }
//---------------------------------------
  public abstract class modeMenuHandler implements ActionListener {
  }
//---------------------------------------
  public abstract class fixedPointMenuHandler implements ActionListener {
  }
//---------------------------------------
  public abstract class thresholdMenuHandler implements ActionListener {
  }
//---------------------------------------
  public abstract class helpMenuHandler implements ActionListener {
  }
//---------------------------------------
  public abstract class planeColChooser implements ActionListener {
  }

} // class SectionViewer
