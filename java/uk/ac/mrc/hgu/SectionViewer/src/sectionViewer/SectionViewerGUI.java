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

/**
 *   GUI for the SectionViewer.
 */
public class SectionViewerGUI extends JPanel {

  /** Toggles output of debugging messages. */
  private final boolean _debug = false;

  /**
   *   Top level container for SectionViewer GUI.
   *   <br>Contains permanentPanel and transientPanel.
   */
  protected JPanel _contentPane;

  // these need to be visible outside of guiInit()

  /**   The preferred height of _contentPane. */
  int totalH = 500;

  /**   The preferred width of _contentPane. */
  int totalW = 400;

  /**   A spacing parameter for elements of the GUI. */
  int bord = 2;

  /**   A spacing parameter for elements of the GUI. */
  int pad = 1;

  /**   A spacing parameter for elements of the GUI. */
  int hgap = 2;

  /**   A spacing parameter for elements of the GUI. */
  int vgap = 2;

  //int titleH;
  /**
   *   Height of the panel containing
   *   <em>pitch yaw roll</em> sliders.
   */
  int pyrH;

  /**
   *   Height of the panel containing <em>fixed line</em> sliders.
   */
  int rotH;

  /**
   *   Height of the transient panel.
   *   <br>The transient panel contains the
   *   <em>pitch yaw roll</em> and <em>fixed line</em> sliders.
   */
  int transH;

  /**   Font used for the menus. */
  //protected final Font _menuFont = new Font("default", Font.PLAIN, 11);
  protected Font _menuFont = null;

  /**   Font used for feedback text. */
  protected final Font feedbackFont = new Font("default", Font.PLAIN, 11);

  /**   Font used for <em>mouse-click anatomy</em> feedback text. */
  protected final Font anatomyFont = new Font("default", Font.ITALIC, 11);

  /**
   *   Outer container for menus.
   *   <br>titleMenuPanel used to contain a text field
   *   but this has been removed.
   */
  JPanel titleMenuPanel = new JPanel();

  /**   Middle container for menus. */
  JPanel menuPanel = new JPanel();

  /**   Inner container for menus. */
  JPanel menuPanel_2 = new JPanel();

  /**   MenuBar for SectionViewer menus. */
  public JMenuBar _menubar = new JMenuBar();

  /**
   *   Container for _imageScrollPane.
   *   <br>Added to feedbackImagePanel.
   */
  JPanel imagePanel = new JPanel();

  /**
   *   Container for instance of WlzImgView.
   *   <br>Added to the Viewport of _imageScrollPane.
   *   <br>(The WlzImgView is added in SectionViewer.)
   */
  JPanel _bigPanel = new JPanel();

  /**
   *   Scrollable container for _bigPanel.
   *   <br>Added to imagePanel.
   */
  protected JScrollPane _imageScrollPane = new JScrollPane();

  /**
   *   Container for imagePanel.
   *   <br>Added to permanentPanel.
   */
  JPanel feedbackImagePanel = new JPanel();

  /**
   *    Java component for setting <em>% magnification</em>
   *    of the image displayed on screen.
   */
  Zoom zoomSetter = new Zoom();

  /**
   *   Java component for setting <em>distance</em> from
   *   the <em>fixed point</em>.
   */
  WSetter distSetter = new WSetter();

  /**
   *   Java component for setting <em>pitch</em> angle.
   */
  WSetter pitchSetter = new WSetter();

  /**
   *   Java component for setting <em>yaw</em> angle.
   */
  WSetter yawSetter = new WSetter();

  /**
   *   Java component for setting <em>roll</em> angle.
   */
  WSetter rollSetter = new WSetter();

  /**
   *   Java component for setting <em>fixed line</em> rotation angle.
   */
  WSetter rotSetter = new WSetter();

  /**
   *   Control for inverting the values of the grey-level image on screen.
   *   <br>Added to invertPanel.
   */
  JButton invertButton = null;

  /**
   *   Container for <em>zoomSetter</em> component, spacerPanel_R
   *   & invertPanel.
   *   <br>Added to basicControlPanel.
   */
  JPanel zoomControlPanel = new JPanel();

  /**  Spacer for <em>zoomSetter</em> and <em>invert</em> controls.  */
  JPanel spacerPanel_R = new JPanel();

  /**
   *   Container for invertButton.
   *   <br>Added to zoomControlPanel.
   */
  JPanel invertPanel = new JPanel();

  /**
   *   Container for zoomControlPanel and <em>distSetter</em>.
   *   <br>Added to permanentPanel.
   */
  JPanel basicControlPanel = new JPanel();

  /**
   *   Container for pitchyawControlPanel & rollControlPanel.
   *   <br>Added to / removed from transientPanel by SectionViewer.
   */
  JPanel pitchYawRollPanel = new JPanel();

  /**
   *   Container for <em>pitchSetter</em> & <em>yawSetter</em>.
   *   <br>Added to pitchYawRollPanel.
   */
  JPanel pitchyawControlPanel = new JPanel();

  /**
   *   Container for <em>rollSetter</em>.
   *   <br>Added to pitchYawRollPanel.
   */
  JPanel rollControlPanel = new JPanel();

  /**
   *   Container for <em>rotSetter</em>.
   *   <br>Added to / removed from transientPanel by SectionViewer.
   */
  JPanel userDefinedRotPanel = new JPanel();

  /*
  JPanel fixedPointControlPanel = new JPanel();
  JPanel fixedPointUserRotPanel = new JPanel();
  */

  //...................

  /**
   *   Container for feedbackImagePanel and basicControlPanel.
   *   <br>Added to _contentPane
   */
  JPanel permanentPanel = new JPanel();

  /**
   *   Container for  pitchYawRollPanel & userDefinedRotPanel.
   *   <br>Added to / removed from _contentPane by SectionViewer.
   */
  JPanel transientPanel = new JPanel();

  //-------------------------------
  /**   file menu name */
  String fileMenuStr = "File";
  /**   file sub-menu name */
  String fileMenu_1str = "Save image";
  /**   file sub-menu name */
  String fileMenu_2str = "Save view settings";
  /**   file sub-menu name */
  String fileMenu_3str = "Load view settings from xml";
  /**   file sub-menu name */
  String fileMenu_4str = "Load view settings from bib";
  /**   file sub-menu name */
  String fileMenu_5str = "Close view";

  /**   file menu */
  JMenu fileMenu = new JMenu(fileMenuStr);

  /**   file sub-menu */
  JMenuItem fileMenu_1 = new JMenuItem(fileMenu_1str);
  /**   file sub-menu */
  JMenuItem fileMenu_2 = new JMenuItem(fileMenu_2str);
  /**   file sub-menu */
  JMenuItem fileMenu_3 = new JMenuItem(fileMenu_3str);
  /**   file sub-menu */
  JMenuItem fileMenu_4 = new JMenuItem(fileMenu_4str);
  /**   file sub-menu */
  JMenuItem fileMenu_5 = new JMenuItem(fileMenu_5str);
  //-------------------------------
  /**   control menu name */
  String controlMenuStr = "Control";

  /**   control menu */
  JMenu controlMenu = new JMenu(controlMenuStr);
  //...................
  /**   control sub-menu name */
  String controlMenu_1str = "Rotation";
  /**   control sub-menu name */
  String controlMenu_1_1str = "Yaw pitch roll";
  /**   control sub-menu name */
  String controlMenu_1_2str = "Fixed_line";
  /**   control sub-menu name */
  String controlMenu_1_3str = "All";

  /**   control sub-menu */
  JMenu controlMenu_1 = new JMenu(controlMenu_1str);
  /**   true if checkBox selected */
  boolean controlMenu_1_1state = false;
  /**   control sub-sub-menu */
  JCheckBoxMenuItem controlMenu_1_1 = new
               JCheckBoxMenuItem(controlMenu_1_1str,
	                         controlMenu_1_1state);
  /**   true if checkBox selected */
  boolean controlMenu_1_2state = false;
  /**   control sub-sub-menu */
  JCheckBoxMenuItem controlMenu_1_2 = new
               JCheckBoxMenuItem(controlMenu_1_2str,
	                         controlMenu_1_2state);
  /**   true if checkBox selected */
  boolean controlMenu_1_3state = false;
  /**   control sub-sub-menu */
  JCheckBoxMenuItem controlMenu_1_3 = new
               JCheckBoxMenuItem(controlMenu_1_3str,
	                         controlMenu_1_3state);
  //...................
  /**   control sub-menu name */
  String controlMenu_2str = "View mode";
  /**   control sub-sub-menu name */
  String controlMenu_2_1str = "Up is up";
  /**   control sub-sub-menu name */
  String controlMenu_2_2str = "Absolute";
  /**   control sub-sub-menu name */
  String controlMenu_2_3str = "Fixed line";

  /**   control sub-menu */
  JMenu controlMenu_2 = new JMenu(controlMenu_2str);
  /**
   *   Mutually exclusive radioButton group.
   *   <br>Contains controlMenu_2_1 controlMenu_2_2 controlMenu_2_3.
   */
  ButtonGroup viewModeGroup = new ButtonGroup();
  /**   control sub-sub-menu */
  JRadioButtonMenuItem controlMenu_2_1 = new
               JRadioButtonMenuItem(controlMenu_2_1str);
  /**   control sub-sub-menu */
  JRadioButtonMenuItem controlMenu_2_2 = new
               JRadioButtonMenuItem(controlMenu_2_2str);
  /**   control sub-sub-menu */
  JRadioButtonMenuItem controlMenu_2_3 = new
               JRadioButtonMenuItem(controlMenu_2_3str);
  //...................
  /**   control sub-menu name */
  String controlMenu_3str = "Fixed point";
  /**   control sub-sub-menu name */
  String controlMenu_3_1str = "Change fixed point using mouse button";
  /**   control sub-sub-menu name */
  String controlMenu_3_2str = "Change fixed point by entering coordinates";
  /**   control sub-sub-menu name */
  String controlMenu_3_3str = "Reset fixed point";

  /**   control sub-menu */
  JMenu controlMenu_3 = new JMenu(controlMenu_3str);
  /**   control sub-sub-menu */
  JMenuItem controlMenu_3_1 = new JMenuItem(controlMenu_3_1str);
  /**   control sub-sub-menu */
  JMenuItem controlMenu_3_2 = new JMenuItem(controlMenu_3_2str);
  /**   control sub-sub-menu */
  JMenuItem controlMenu_3_3 = new JMenuItem(controlMenu_3_3str);
  //...................
  /**   control sub-menu name */
  String controlMenu_4str = "Fixed line end-point";
  /**   control sub-sub-menu name */
  String controlMenu_4_1str = "Change fixed line using mouse button";
  /**   control sub-sub-menu name */
  String controlMenu_4_2str = "Change fixed line by entering coordinates";
  /**   control sub-sub-menu name */
  String controlMenu_4_3str = "Reset fixed line";

  /**   control sub-menu */
  JMenu controlMenu_4 = new JMenu(controlMenu_4str);
  /**   control sub-sub-menu */
  JMenuItem controlMenu_4_1 = new JMenuItem(controlMenu_4_1str);
  /**   control sub-sub-menu */
  JMenuItem controlMenu_4_2 = new JMenuItem(controlMenu_4_2str);
  /**   control sub-sub-menu */
  JMenuItem controlMenu_4_3 = new JMenuItem(controlMenu_4_3str);
  //...................
  /**   control sub-menu name */
  String controlMenu_5str = "Reset controls";
  /**   control sub-menu */
  JMenuItem controlMenu_5 = new JMenuItem(controlMenu_5str);
  //-------------------------------
  /**   show menu name */
  String showMenuStr = "Show";
  /**   show menu */
  JMenu showMenu = new JMenu(showMenuStr);

  /**   show sub-menu name */
  String showMenu_1str = "Cursor feedback";
  /**   true if checkBox selected */
  boolean showMenu_1state = true;
  /**   show sub-menu */
  JCheckBoxMenuItem showMenu_1 = new
               JCheckBoxMenuItem(showMenu_1str,
                                 showMenu_1state);

  /**   show sub-menu name */
  String showMenu_2str = "Intersection of views";
  /**   true if checkbox selected */
  boolean showMenu_2state = false;
  /**   show sub-menu */
  JCheckBoxMenuItem showMenu_2 = new
               JCheckBoxMenuItem(showMenu_2str,
                                 showMenu_2state);

  /**   show sub-menu name */
  String showMenu_3str = "Mouse-click anatomy";
  /**   true if checkbox selected */
  boolean showMenu_3state = false;
  /**   show sub-menu */
  JCheckBoxMenuItem showMenu_3 = new
               JCheckBoxMenuItem(showMenu_3str,
                                 showMenu_3state);

  /**   show sub-menu name */
  String showMenu_4str = "Fixed point";
  /**   true if checkbox selected */
  boolean showMenu_4state = false;
  /**   show sub-menu */
  JCheckBoxMenuItem showMenu_4 = new
               JCheckBoxMenuItem(showMenu_4str,
                                 showMenu_4state);

  /**   show sub-menu name */
  String showMenu_5str = "Fixed line";
  /**   true if checkbox selected */
  boolean showMenu_5state = false;
  /**   show sub-menu */
  JCheckBoxMenuItem showMenu_5 = new
               JCheckBoxMenuItem(showMenu_5str,
                                 showMenu_5state);
  //-------------------------------
  // not used at present, but don't delete
  /**   threshold menu name */
  String thresholdMenuStr = "Threshold";
  /**   threshold sub-menu name */
  String thresholdMenu_1str = "Enable constraint definition";
  /**   threshold sub-menu name */
  String thresholdMenu_2str = "Enable thresholding";
  /**   threshold sub-menu name */
  String thresholdMenu_3str = "Show threshold constraint";
  /**   threshold sub-menu name */
  String thresholdMenu_4str = "Remove threshold constraint";

  /**   true if checkbox selected */
  boolean thresholdMenu_1state = false;
  /**   true if checkbox selected */
  boolean thresholdMenu_2state = false;
  /**   true if checkbox selected */
  boolean thresholdMenu_3state = false;

  /**   threshold menu */
  JMenu thresholdMenu = new JMenu(thresholdMenuStr);
  /**   threshold sub-menu */
  JCheckBoxMenuItem thresholdMenu_1 = new
               JCheckBoxMenuItem(thresholdMenu_1str, thresholdMenu_1state);
  /**   threshold sub-menu */
  JCheckBoxMenuItem thresholdMenu_2 = new
               JCheckBoxMenuItem(thresholdMenu_2str, thresholdMenu_2state);
  /**   threshold sub-menu */
  JCheckBoxMenuItem thresholdMenu_3 = new
               JCheckBoxMenuItem(thresholdMenu_3str, thresholdMenu_3state);
  /**   threshold sub-menu */
  JMenuItem thresholdMenu_4 = new JMenuItem(thresholdMenu_4str);

  //-------------------------------
  /**   help menu name */
  String helpMenuStr = "Help";
  /**   help sub-menu name */
  String helpMenu_1str = "Contents";
  /**   help sub-menu name */
  String helpMenu_2str = "Index";
  /**   help sub-menu name */
  String helpMenu_3str = "Search";
  /**   help sub-menu */
  JMenuItem helpMenu_1 = new JMenuItem(helpMenu_1str);
  /**   help sub-menu */
  JMenuItem helpMenu_2 = new JMenuItem(helpMenu_2str);
  /**   help sub-menu */
  JMenuItem helpMenu_3 = new JMenuItem(helpMenu_3str);

  /**   help menu */
  public JMenu helpMenu = new JMenu(helpMenuStr);
  //-------------------------------
  //-------------------------------
  /**
   *   Container for feedbackPanel_3, fbAnatPanel & secColorClt.
   *   <br>Added to / Removed from feedbackImagePanel by SectionViewer.
   */
  JPanel feedbackPanel = new JPanel();

  /**
   *   Container for xyzTextField.
   *   <br>Added to feedbackPanel_3.
   */
  JPanel feedbackPanel_1 = new JPanel();

  /**
   *   Spacer between <em>positional</em>
   *   & <em>grey value</em> feedback.
   *   <br>Added to feedbackPanel_2.
   */
  JPanel feedbackPanel_1a = new JPanel();

  /**
   *   Container for feedbackPanel_1a & valueTextField.
   *   <br>Added to feedbackPanel_3.
   */
  JPanel feedbackPanel_2 = new JPanel();

  /**
   *   Container for feedbackPanel_1 & feedbackPanel_2.
   *   <br>Added to feedbackPanel.
   */
  JPanel feedbackPanel_3 = new JPanel();

  /*  Unused Containers.
  JPanel feedbackPanel_4 = new JPanel();
  JPanel feedbackPanel_5 = new JPanel();
  JPanel feedbackPanel_6 = new JPanel();
  JPanel feedbackPanel_7 = new JPanel();
  */

  /**
   *   Container for anatomyTextField.
   *   <br>Added to feedbackPanel.
   *   <br>Must to be added to BorderLayout.CENTER to allow for expansion.
   */
  JPanel fbAnatPanel = new JPanel();

  /**
   *   Display for <em>positional (x,y,z)</em> feedback text.
   */
  protected JTextField xyzTextField = new JTextField();

  /**   Display for <em>grey value</em> at cursor. */
  JTextField valueTextField = new JTextField();

  /**   Display for name of <em>mouse-click anatomy</em> at cursor. */
  ScrollableTextField anatomyTextField = new ScrollableTextField();

  /* unused variables
  String titleString = "Ref:";
  Color titleCol = new Color(200, 220, 220);
  */

  /**   Background colour for feedback text fields. */
  Color fbCol = new Color(220, 220, 200);

  /**
   *   Initiates colour chooser dialogue for SectionViewer.
   *   <br>Added to feedbackPanel.
   */
  public JButton secColorClt = new JButton("");

  /**
   *   Indicator for the number of open SectionViewers.
   *   <br>Used to assign a colour to a newly opened SectionViewer.
   */
  public static int nSV = 0;

  /**
   *   List of colours to be assigned initially to SectionViewers.
   *   The list is repeated if required.
   */
  private Color[] planeColor = new Color[] {Color.red, Color.yellow,
      Color.blue, Color.pink, Color.green, Color.cyan};

  //=========================================================
  // constructor
  //=========================================================
  /**
   *   Constructs the GUI for a SectionViewer.
   */
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
  /**
   *   Most of the work of constucting the GUI is done here..
   */
  private void guiInit() throws Exception {

    if(_debug) System.out.println("enter guiInit");

    _contentPane = this;

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
    int invertH = 20;
    int invertW = 20;

    int permH;
    //......................................
    int pyH = 40;
    int rH = 20;
    pyrH = pyH+rH;
    rotH = 25;
    transH = pyrH+vgap+rotH;

//----------------------------------------------------------
    /*
    fileMenu.setFont(_menuFont);
    fileMenu_1.setFont(_menuFont);
    fileMenu_2.setFont(_menuFont);
    fileMenu_3.setFont(_menuFont);
    fileMenu_4.setFont(_menuFont);
    fileMenu_5.setFont(_menuFont);
    */
    _menuFont = _menubar.getFont();

    fileMenu.add(fileMenu_1);
    fileMenu.addSeparator();
    fileMenu.add(fileMenu_2);
    fileMenu.add(fileMenu_3);
    fileMenu.add(fileMenu_4);
    fileMenu.addSeparator();
    fileMenu.add(fileMenu_5);
    fileMenu_1.setEnabled(true);
    fileMenu_2.setEnabled(true);
    fileMenu_3.setEnabled(true);
    fileMenu_4.setEnabled(true);
    fileMenu_5.setEnabled(true);
 //...............................
 //...............................
    /*
    controlMenu.setFont(_menuFont);
    controlMenu_1.setFont(_menuFont);
    controlMenu_2.setFont(_menuFont);
    controlMenu_3.setFont(_menuFont);
    controlMenu_4.setFont(_menuFont);
    controlMenu_5.setFont(_menuFont);
    controlMenu_1_1.setFont(_menuFont);
    controlMenu_1_2.setFont(_menuFont);
    controlMenu_1_3.setFont(_menuFont);
    controlMenu_2_1.setFont(_menuFont);
    controlMenu_2_2.setFont(_menuFont);
    controlMenu_2_3.setFont(_menuFont);
    controlMenu_3_1.setFont(_menuFont);
    controlMenu_3_2.setFont(_menuFont);
    controlMenu_3_3.setFont(_menuFont);
    controlMenu_4_1.setFont(_menuFont);
    controlMenu_4_2.setFont(_menuFont);
    controlMenu_4_3.setFont(_menuFont);
    */

    controlMenu.add(controlMenu_1);
    controlMenu.add(controlMenu_2);
    controlMenu.add(controlMenu_3);
    controlMenu.add(controlMenu_4);
    controlMenu.add(controlMenu_5);

    // rotation
    controlMenu_1.add(controlMenu_1_1);
    controlMenu_1.add(controlMenu_1_2);
    controlMenu_1.add(controlMenu_1_3);

    // view mode
    controlMenu_2.add(controlMenu_2_1);
    controlMenu_2.add(controlMenu_2_2);
    controlMenu_2.add(controlMenu_2_3);
    viewModeGroup.add(controlMenu_2_1);
    viewModeGroup.add(controlMenu_2_2);
    viewModeGroup.add(controlMenu_2_3);

    controlMenu_2_1.setSelected(true);
    controlMenu_2_3.setEnabled(false);

    // fixed point
    controlMenu_3.add(controlMenu_3_1);
    controlMenu_3.add(controlMenu_3_2);
    controlMenu_3.add(controlMenu_3_3);

    // fixed line
    controlMenu_4.add(controlMenu_4_1);
    controlMenu_4.add(controlMenu_4_2);
    controlMenu_4.add(controlMenu_4_3);
 //...............................
 //...............................
    /*
    showMenu.setFont(_menuFont);
    showMenu_1.setFont(_menuFont);
    showMenu_2.setFont(_menuFont);
    showMenu_3.setFont(_menuFont);
    showMenu_4.setFont(_menuFont);
    showMenu_5.setFont(_menuFont);
    */

    showMenu.add(showMenu_1); // cursor feedback
    showMenu.add(showMenu_2); // intersection of views
    showMenu.add(showMenu_3); // mouse-click anatomy
    showMenu.add(showMenu_4); // fixed point
    showMenu.add(showMenu_5); // fixed line
    //showMenu_5.setEnabled(false);
 //...............................
 //...............................
    /*
    thresholdMenu.setFont(_menuFont);
    thresholdMenu_1.setFont(_menuFont);
    thresholdMenu_2.setFont(_menuFont);
    thresholdMenu_3.setFont(_menuFont);
    thresholdMenu_4.setFont(_menuFont);
    */

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
    /*
    helpMenu.setFont(_menuFont);
    helpMenu_1.setFont(_menuFont);
    helpMenu_2.setFont(_menuFont);
    helpMenu_3.setFont(_menuFont);
    */

    helpMenu.add(helpMenu_1);
    helpMenu.add(helpMenu_2);
    helpMenu.add(helpMenu_3);
    helpMenu_1.setEnabled(true);
    helpMenu_2.setEnabled(true);
    helpMenu_3.setEnabled(true);
    //----------------------------------------------------------
    _menubar.add(fileMenu);
    _menubar.add(controlMenu);
    _menubar.add(showMenu);
    //_menubar.add(thresholdMenu);
    _menubar.add(helpMenu);
    //----------------------------------------------------------
    menuH = _menubar.getFontMetrics(_menuFont).getHeight() + bord;
    //......................................
    fbH = xyzTextField.getFontMetrics(feedbackFont).getHeight() + bord;
    //......................................
    xyzH = fbH+2*bord;
    valH = fbH+2*bord;
    anatH = fbH+2*bord;
    fbimgH = fbH+imgH+vgap;
    permH = menuH+vgap+fbimgH+vgap+distH;
    //----------------------------------------------------------

    xyzTextField.setFont(feedbackFont);
    xyzTextField.setEditable(false);
    xyzTextField.setBackground(fbCol);

    valueTextField.setFont(feedbackFont);
    valueTextField.setEditable(false);
    valueTextField.setBackground(fbCol);

    anatomyTextField.setFont(feedbackFont);
    anatomyTextField.setEditable(false);
    anatomyTextField.setBackground(fbCol);

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
    //secColorClt.setPreferredSize(new Dimension(titleH, titleH));
    secColorClt.setPreferredSize(new Dimension(fbH, fbH));
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

    invertButton = new JButton();
    invertButton.setPreferredSize(new Dimension(invertW, invertH));
    invertButton.setBackground(Color.white);
    invertButton.setForeground(Color.black);
    invertPanel.setLayout(new BorderLayout(hgap, vgap));
    invertPanel.add(invertButton, BorderLayout.EAST);

    distSetter.setBgc(new Color(255, 255, 230));
    distSetter.setLabelWidth(60);
    distSetter.setTextWidth(40);
    distSetter.setSliderLabel("dist");

    zoomControlPanel.setLayout(new BorderLayout(hgap, vgap));
    zoomControlPanel.add(zoomSetter, BorderLayout.WEST);
    zoomControlPanel.add(spacerPanel_R, BorderLayout.CENTER);
    zoomControlPanel.add(invertPanel, BorderLayout.EAST);

    basicControlPanel.setPreferredSize(new Dimension(totalW, 2*distH));
    basicControlPanel.setLayout(new BorderLayout(hgap, vgap));
    basicControlPanel.setBorder(BorderFactory.createCompoundBorder(
                 BorderFactory.createEmptyBorder(3,1,3,1),
		             BorderFactory.createEtchedBorder(EtchedBorder.LOWERED)));
    basicControlPanel.add(zoomControlPanel, BorderLayout.NORTH);
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
    yawSetter.setMax(360.0);
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
    rollSetter.setMax(360.0);
    rollSetter.setValue(0.0);

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
    rotSetter.setSliderLabel("fixed line");
    rotSetter.setMin(-180.0);
    rotSetter.setMax(180.0);
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

    _contentPane.setLayout(new BorderLayout());
    _contentPane.setPreferredSize(new Dimension(totalW+pad, totalH+pad));
    _contentPane.setMinimumSize(new Dimension(totalW+pad, totalH+pad));
    _contentPane.add(permanentPanel, BorderLayout.CENTER);
    _contentPane.add(transientPanel, BorderLayout.SOUTH);

    //----------------------------------------------------------
    if(_debug) System.out.println("exit guiInit");
  } // guiInit()


//----------------------------------------------------------

//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class fileMenuHandler implements ActionListener {
  }
//---------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class controlMenuHandler implements ActionListener {
  }
//---------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class modeMenuHandler implements ActionListener {
  }
//---------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class fixedPointMenuHandler implements ActionListener {
  }
//---------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class thresholdMenuHandler implements ActionListener {
  }
//---------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class helpMenuHandler implements ActionListener {
  }
//---------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class planeColChooser implements ActionListener {
  }
//---------------------------------------
  /**  Declaration for fileMenuHandler. */
  public abstract class invertButtonHandler implements ActionListener {
  }

} // class SectionViewerGUI
