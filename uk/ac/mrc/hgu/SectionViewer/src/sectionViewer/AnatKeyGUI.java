package sectionViewer;
import sectionViewer.*;

import java.util.*;
import java.net.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;

public class AnatKeyGUI extends JFrame{

   private String SLASH = System.getProperty("file.separator");

   static public int _cols[];
   static protected int _nrows = 6;
   /*
      colours are taken from:
      src/Applications/MAPaint/MAPaintResources.h
      to be consistent with Editorial Office
      (note however that red has been moved from #6 to #2
      and is solid for E.O.)
   */
   static protected int _rgbt[][] = {{128,96,0,125},
                                     {255,0,0,125},
                                     {255,255,0,125},
                                     {0,0,255,125},
                                     {0,255,0,125},
                                     {0,255,255,125}};

   private JPanel topPanel = null;

   static protected JButton colBtn0 = null;
   static protected JButton colBtn1 = null;
   static protected JButton colBtn2 = null;
   static protected JButton colBtn3 = null;
   static protected JButton colBtn4 = null;
   static protected JButton colBtn5 = null;

   static protected JButton btn00 = null;
   static protected JButton btn10 = null;
   static protected JButton btn20 = null;
   static protected JButton btn30 = null;
   static protected JButton btn40 = null;
   static protected JButton btn50 = null;

   static protected JButton btn01 = null;
   static protected JButton btn11 = null;
   static protected JButton btn21 = null;
   static protected JButton btn31 = null;
   static protected JButton btn41 = null;
   static protected JButton btn51 = null;

   static protected ScrollableTextField tf0 = null;
   static protected ScrollableTextField tf1 = null;
   static protected ScrollableTextField tf2 = null;
   static protected ScrollableTextField tf3 = null;
   static protected ScrollableTextField tf4 = null;
   static protected ScrollableTextField tf5 = null;

   static protected JCheckBox cbx0 = null;
   static protected JCheckBox cbx1 = null;
   static protected JCheckBox cbx2 = null;
   static protected JCheckBox cbx3 = null;
   static protected JCheckBox cbx4 = null;
   static protected JCheckBox cbx5 = null;

   static protected boolean cbx0State;
   static protected boolean cbx1State;
   static protected boolean cbx2State;
   static protected boolean cbx3State;
   static protected boolean cbx4State;
   static protected boolean cbx5State;

   static protected JButton btn02 = null;
   static protected JButton btn12 = null;
   static protected JButton btn22 = null;
   static protected JButton btn32 = null;
   static protected JButton btn42 = null;
   static protected JButton btn52 = null;

   static protected ImageIcon leftIcon = null;
   static protected ImageIcon rightIcon = null;
   static protected ImageIcon zapIcon = null;
   static protected ImageIcon replaceIcon = null;

   static protected JToolTip colTip0 = null;
   static protected JToolTip btnTip00 = null;
   static protected JToolTip btnTip01 = null;
   static protected JToolTip btnTip02 = null;
   static protected JToolTip cbxTip0 = null;

   static protected String colTipText = "choose colour";
   static protected String leftTipText = "start of text";
   static protected String rightTipText = "end of text";
   //static protected String textFieldTipText = "drag mouse to scroll text";
   static protected String visTipText = "toggle visibility";
   static protected String zapTipText = "remove from key";

   static Color vis;
   static Color notvis;
//-------------------------------------------------------------
   protected AnatKeyGUI(String str) {
      super(str);
      setCols();
      topPanel = (JPanel)this.getContentPane();
      topPanel.setLayout(new BorderLayout());
      topPanel.setPreferredSize(new Dimension(250,150));
      cbx0State = false;
      cbx1State = false;
      cbx2State = false;
      cbx3State = false;
      cbx4State = false;
      cbx5State = false;
      makeGUI();
   }

//-------------------------------------------------------------
   protected void setCols() {

      // make initial color scheme consistent with MAPaint etc

      _cols = new int[_nrows];
      _cols[0] = (_rgbt[0][2])|
                 (_rgbt[0][1]<<8)|
                 (_rgbt[0][0]<<16)|
                 (_rgbt[0][3]<<24);
      _cols[1] = (_rgbt[1][2])|
                 (_rgbt[1][1]<<8)|
                 (_rgbt[1][0]<<16)|
                 (_rgbt[1][3]<<24);
      _cols[2] = (_rgbt[2][2])|
                 (_rgbt[2][1]<<8)|
                 (_rgbt[2][0]<<16)|
                 (_rgbt[2][3]<<24);
      _cols[3] = (_rgbt[3][2])|
                 (_rgbt[3][1]<<8)|
                 (_rgbt[3][0]<<16)|
                 (_rgbt[3][3]<<24);
      _cols[4] = (_rgbt[4][2])|
                 (_rgbt[4][1]<<8)|
                 (_rgbt[4][0]<<16)|
                 (_rgbt[4][3]<<24);
      _cols[5] = (_rgbt[5][2])|
                 (_rgbt[5][1]<<8)|
                 (_rgbt[5][0]<<16)|
                 (_rgbt[5][3]<<24);
   }
//-------------------------------------------------------------
   private void makeGUI() {

      boolean hasAlpha = false;

      int colW = 20;
      int btnW1 = 12;
      int btnW2 = 20;
      int btnW3 = 25;
      int cbW = 20;
      int gap = 2;
      int gap1 = 1;

      String imgPath = "sectionViewer"+SLASH+"images"+SLASH;
      leftIcon = new ImageIcon(imgPath+"left.gif");
      rightIcon = new ImageIcon(imgPath+"right.gif");
      zapIcon = new ImageIcon(imgPath+"zap.gif");
      replaceIcon = new ImageIcon(imgPath+"replace.gif");

      vis = new Color(250,250,250);
      notvis = new Color(200,200,200);

      JPanel kTopPanel = new JPanel();
      kTopPanel.setLayout(new GridLayout(_nrows,1));

//......................................
      JPanel kPanel0 = new JPanel(); // 1 row of the key
      JPanel kPanel00 = new JPanel(); // holds colour panel & arrow button
      JPanel kPanel01 = new JPanel(); // holds chkbx and remove button
      kPanel0.setLayout(new BorderLayout(gap,gap));
      kPanel00.setLayout(new BoxLayout(kPanel00, BoxLayout.X_AXIS));
      kPanel01.setLayout(new BorderLayout(gap1,gap1));
      kPanel00.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      kPanel01.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      colBtn0 = new JButton();
      colBtn0.setBackground(new Color(_cols[0], hasAlpha));
      JPanel colPanel0 = new JPanel(new GridLayout());
      colPanel0.add(colBtn0);
      colPanel0.setPreferredSize(new Dimension(colW,0));
      JPanel txtPanel0 = new JPanel();
      txtPanel0.setLayout(new BorderLayout(gap1,gap1));
      txtPanel0.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      tf0 = new ScrollableTextField("", 50);
      tf0.setEditable(false);
      tf0.setBackground(notvis);
      txtPanel0.add(tf0, BorderLayout.CENTER);
      JPanel navPanel0 = new JPanel();
      navPanel0.setPreferredSize(new Dimension(btnW3,0));
      navPanel0.setLayout(new BoxLayout(navPanel0, BoxLayout.X_AXIS));
      JPanel btnPanel00 = new JPanel();
      btnPanel00.setPreferredSize(new Dimension(btnW1,0));
      btnPanel00.setLayout(new BorderLayout(gap1,gap1));
      btn00 = new JButton();
      if(leftIcon != null) {
	 btn00.setIcon(leftIcon);
      }
      btnPanel00.add(btn00, BorderLayout.CENTER);
      JPanel btnPanel01 = new JPanel();
      btnPanel01.setPreferredSize(new Dimension(btnW1,0));
      btnPanel01.setLayout(new BorderLayout(gap1,gap1));
      btn01 = new JButton();
      if(leftIcon != null) {
	 btn01.setIcon(rightIcon);
      }
      btnPanel01.add(btn01, BorderLayout.CENTER);
      navPanel0.add(btnPanel00);
      navPanel0.add(btnPanel01);

      JPanel cbxPanel0 = new JPanel();
      cbxPanel0.setPreferredSize(new Dimension(btnW1,0));
      cbxPanel0.setLayout(new FlowLayout(FlowLayout.CENTER, gap1,gap1));
      cbx0 = new JCheckBox();
      cbxPanel0.add(cbx0, BorderLayout.CENTER);
      cbx0.setSelected(cbx0State);
      JPanel btnPanel02 = new JPanel();
      btnPanel02.setPreferredSize(new Dimension(btnW2,0));
      btnPanel02.setLayout(new BorderLayout(gap1,gap1));
      btn02 = new JButton();
      btnPanel02.add(btn02, BorderLayout.CENTER);

      kPanel00.add(colPanel0);
      kPanel00.add(navPanel0);
      kPanel01.add(cbxPanel0, BorderLayout.CENTER);
      kPanel01.add(btnPanel02, BorderLayout.EAST);
      kPanel0.add(kPanel00, BorderLayout.WEST);
      kPanel0.add(txtPanel0, BorderLayout.CENTER);
      kPanel0.add(kPanel01, BorderLayout.EAST);

//......................................
      JPanel kPanel1 = new JPanel(); // 1 row of the key
      JPanel kPanel10 = new JPanel(); // holds colour panel & arrow button
      JPanel kPanel11 = new JPanel(); // holds chkbx and remove button
      kPanel1.setLayout(new BorderLayout(gap,gap));
      kPanel10.setLayout(new BoxLayout(kPanel10, BoxLayout.X_AXIS));
      kPanel11.setLayout(new BorderLayout(gap1,gap1));
      kPanel10.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      kPanel11.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      colBtn1 = new JButton();
      colBtn1.setBackground(new Color(_cols[1], hasAlpha));
      JPanel colPanel1 = new JPanel(new GridLayout());
      colPanel1.add(colBtn1);
      colPanel1.setPreferredSize(new Dimension(colW,0));
      JPanel txtPanel1 = new JPanel();
      txtPanel1.setLayout(new BorderLayout(gap1,gap1));
      txtPanel1.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      tf1 = new ScrollableTextField("", 50);
      tf1.setEditable(false);
      tf1.setBackground(notvis);
      txtPanel1.add(tf1, BorderLayout.CENTER);
      JPanel navPanel1 = new JPanel();
      navPanel1.setPreferredSize(new Dimension(btnW3,0));
      navPanel1.setLayout(new BoxLayout(navPanel1, BoxLayout.X_AXIS));
      JPanel btnPanel10 = new JPanel();
      btnPanel10.setPreferredSize(new Dimension(btnW1,0));
      btnPanel10.setLayout(new BorderLayout(gap1,gap1));
      btn10 = new JButton();
      if(leftIcon != null) {
	 btn10.setIcon(leftIcon);
      }
      btnPanel10.add(btn10, BorderLayout.CENTER);
      JPanel btnPanel11 = new JPanel();
      btnPanel11.setPreferredSize(new Dimension(btnW1,0));
      btnPanel11.setLayout(new BorderLayout(gap1,gap1));
      btn11 = new JButton();
      if(leftIcon != null) {
	 btn11.setIcon(rightIcon);
      }
      btnPanel11.add(btn11, BorderLayout.CENTER);
      navPanel1.add(btnPanel10);
      navPanel1.add(btnPanel11);
      JPanel cbxPanel1 = new JPanel();
      cbxPanel1.setPreferredSize(new Dimension(btnW1,0));
      cbxPanel1.setLayout(new FlowLayout(FlowLayout.CENTER, gap1,gap1));
      cbx1 = new JCheckBox();
      cbxPanel1.add(cbx1, BorderLayout.CENTER);
      cbx1.setSelected(cbx1State);
      JPanel btnPanel12 = new JPanel();
      btnPanel12.setPreferredSize(new Dimension(btnW2,0));
      btnPanel12.setLayout(new BorderLayout(gap1,gap1));
      btn12 = new JButton();
      btnPanel12.add(btn12, BorderLayout.CENTER);

      kPanel10.add(colPanel1);
      kPanel10.add(navPanel1);
      kPanel11.add(cbxPanel1, BorderLayout.CENTER);
      kPanel11.add(btnPanel12, BorderLayout.EAST);
      kPanel1.add(kPanel10, BorderLayout.WEST);
      kPanel1.add(txtPanel1, BorderLayout.CENTER);
      kPanel1.add(kPanel11, BorderLayout.EAST);

//......................................
      JPanel kPanel2 = new JPanel(); // 1 row of the key
      JPanel kPanel20 = new JPanel(); // holds colour panel & arrow button
      JPanel kPanel21 = new JPanel(); // holds chkbx and remove button
      kPanel2.setLayout(new BorderLayout(gap,gap));
      kPanel20.setLayout(new BoxLayout(kPanel20, BoxLayout.X_AXIS));
      kPanel21.setLayout(new BorderLayout(gap1,gap1));
      kPanel20.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      kPanel21.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      colBtn2 = new JButton();
      colBtn2.setBackground(new Color(_cols[2], hasAlpha));
      JPanel colPanel2 = new JPanel(new GridLayout());
      colPanel2.add(colBtn2);
      colPanel2.setPreferredSize(new Dimension(colW,0));
      JPanel txtPanel2 = new JPanel();
      txtPanel2.setLayout(new BorderLayout(gap1,gap1));
      txtPanel2.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      tf2 = new ScrollableTextField("", 50);
      tf2.setEditable(false);
      tf2.setBackground(notvis);
      txtPanel2.add(tf2, BorderLayout.CENTER);
      JPanel navPanel2 = new JPanel();
      navPanel2.setPreferredSize(new Dimension(btnW3,0));
      navPanel2.setLayout(new BoxLayout(navPanel2, BoxLayout.X_AXIS));
      JPanel btnPanel20 = new JPanel();
      btnPanel20.setPreferredSize(new Dimension(btnW1,0));
      btnPanel20.setLayout(new BorderLayout(gap1,gap1));
      btn20 = new JButton();
      if(leftIcon != null) {
	 btn20.setIcon(leftIcon);
      }
      btnPanel20.add(btn20, BorderLayout.CENTER);
      JPanel btnPanel21 = new JPanel();
      btnPanel21.setPreferredSize(new Dimension(btnW1,0));
      btnPanel21.setLayout(new BorderLayout(gap1,gap1));
      btn21 = new JButton();
      if(leftIcon != null) {
	 btn21.setIcon(rightIcon);
      }
      btnPanel21.add(btn21, BorderLayout.CENTER);
      navPanel2.add(btnPanel20);
      navPanel2.add(btnPanel21);
      JPanel cbxPanel2 = new JPanel();
      cbxPanel2.setPreferredSize(new Dimension(btnW1,0));
      cbxPanel2.setLayout(new FlowLayout(FlowLayout.CENTER, gap1,gap1));
      cbx2 = new JCheckBox();
      cbxPanel2.add(cbx2, BorderLayout.CENTER);
      cbx2.setSelected(cbx2State);
      JPanel btnPanel22 = new JPanel();
      btnPanel22.setPreferredSize(new Dimension(btnW2,0));
      btnPanel22.setLayout(new BorderLayout(gap1,gap1));
      btn22 = new JButton();
      btnPanel22.add(btn22, BorderLayout.CENTER);

      kPanel20.add(colPanel2);
      kPanel20.add(navPanel2);
      kPanel21.add(cbxPanel2, BorderLayout.CENTER);
      kPanel21.add(btnPanel22, BorderLayout.EAST);
      kPanel2.add(kPanel20, BorderLayout.WEST);
      kPanel2.add(txtPanel2, BorderLayout.CENTER);
      kPanel2.add(kPanel21, BorderLayout.EAST);

//......................................
      JPanel kPanel3 = new JPanel(); // 1 row of the key
      JPanel kPanel30 = new JPanel(); // holds colour panel & arrow button
      JPanel kPanel31 = new JPanel(); // holds chkbx and remove button
      kPanel3.setLayout(new BorderLayout(gap,gap));
      kPanel30.setLayout(new BoxLayout(kPanel30, BoxLayout.X_AXIS));
      kPanel31.setLayout(new BorderLayout(gap1,gap1));
      kPanel30.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      kPanel31.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      colBtn3 = new JButton();
      colBtn3.setBackground(new Color(_cols[3], hasAlpha));
      JPanel colPanel3 = new JPanel(new GridLayout());
      colPanel3.add(colBtn3);
      colPanel3.setPreferredSize(new Dimension(colW,0));
      JPanel txtPanel3 = new JPanel();
      txtPanel3.setLayout(new BorderLayout(gap1,gap1));
      txtPanel3.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      tf3 = new ScrollableTextField("", 50);
      tf3.setEditable(false);
      tf3.setBackground(notvis);
      txtPanel3.add(tf3, BorderLayout.CENTER);
      JPanel navPanel3 = new JPanel();
      navPanel3.setPreferredSize(new Dimension(btnW3,0));
      navPanel3.setLayout(new BoxLayout(navPanel3, BoxLayout.X_AXIS));
      JPanel btnPanel30 = new JPanel();
      btnPanel30.setPreferredSize(new Dimension(btnW1,0));
      btnPanel30.setLayout(new BorderLayout(gap1,gap1));
      btn30 = new JButton();
      if(leftIcon != null) {
	 btn30.setIcon(leftIcon);
      }
      btnPanel30.add(btn30, BorderLayout.CENTER);
      JPanel btnPanel31 = new JPanel();
      btnPanel31.setPreferredSize(new Dimension(btnW1,0));
      btnPanel31.setLayout(new BorderLayout(gap1,gap1));
      btn31 = new JButton();
      if(leftIcon != null) {
	 btn31.setIcon(rightIcon);
      }
      btnPanel31.add(btn31, BorderLayout.CENTER);
      navPanel3.add(btnPanel30);
      navPanel3.add(btnPanel31);
      JPanel cbxPanel3 = new JPanel();
      cbxPanel3.setPreferredSize(new Dimension(btnW1,0));
      cbxPanel3.setLayout(new FlowLayout(FlowLayout.CENTER, gap1,gap1));
      cbx3 = new JCheckBox();
      cbxPanel3.add(cbx3, BorderLayout.CENTER);
      cbx3.setSelected(cbx3State);
      JPanel btnPanel32 = new JPanel();
      btnPanel32.setPreferredSize(new Dimension(btnW2,0));
      btnPanel32.setLayout(new BorderLayout(gap1,gap1));
      btn32 = new JButton();
      btnPanel32.add(btn32, BorderLayout.CENTER);

      kPanel30.add(colPanel3);
      kPanel30.add(navPanel3);
      kPanel31.add(cbxPanel3, BorderLayout.CENTER);
      kPanel31.add(btnPanel32, BorderLayout.EAST);
      kPanel3.add(kPanel30, BorderLayout.WEST);
      kPanel3.add(txtPanel3, BorderLayout.CENTER);
      kPanel3.add(kPanel31, BorderLayout.EAST);

//......................................
      JPanel kPanel4 = new JPanel(); // 1 row of the key
      JPanel kPanel40 = new JPanel(); // holds colour panel & arrow button
      JPanel kPanel41 = new JPanel(); // holds chkbx and remove button
      kPanel4.setLayout(new BorderLayout(gap,gap));
      kPanel40.setLayout(new BoxLayout(kPanel40, BoxLayout.X_AXIS));
      kPanel41.setLayout(new BorderLayout(gap1,gap1));
      kPanel40.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      kPanel41.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      colBtn4 = new JButton();
      colBtn4.setBackground(new Color(_cols[4], hasAlpha));
      JPanel colPanel4 = new JPanel(new GridLayout());
      colPanel4.add(colBtn4);
      colPanel4.setPreferredSize(new Dimension(colW,0));
      JPanel txtPanel4 = new JPanel();
      txtPanel4.setLayout(new BorderLayout(gap1,gap1));
      txtPanel4.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      tf4 = new ScrollableTextField("", 50);
      tf4.setEditable(false);
      tf4.setBackground(notvis);
      txtPanel4.add(tf4, BorderLayout.CENTER);
      JPanel navPanel4 = new JPanel();
      navPanel4.setPreferredSize(new Dimension(btnW3,0));
      navPanel4.setLayout(new BoxLayout(navPanel4, BoxLayout.X_AXIS));
      JPanel btnPanel40 = new JPanel();
      btnPanel40.setPreferredSize(new Dimension(btnW1,0));
      btnPanel40.setLayout(new BorderLayout(gap1,gap1));
      btn40 = new JButton();
      if(leftIcon != null) {
	 btn40.setIcon(leftIcon);
      }
      btnPanel40.add(btn40, BorderLayout.CENTER);
      JPanel btnPanel41 = new JPanel();
      btnPanel41.setPreferredSize(new Dimension(btnW1,0));
      btnPanel41.setLayout(new BorderLayout(gap1,gap1));
      btn41 = new JButton();
      if(leftIcon != null) {
	 btn41.setIcon(rightIcon);
      }
      btnPanel41.add(btn41, BorderLayout.CENTER);
      navPanel4.add(btnPanel40);
      navPanel4.add(btnPanel41);
      JPanel cbxPanel4 = new JPanel();
      cbxPanel4.setPreferredSize(new Dimension(btnW1,0));
      cbxPanel4.setLayout(new FlowLayout(FlowLayout.CENTER, gap1,gap1));
      cbx4 = new JCheckBox();
      cbxPanel4.add(cbx4, BorderLayout.CENTER);
      cbx4.setSelected(cbx4State);
      JPanel btnPanel42 = new JPanel();
      btnPanel42.setPreferredSize(new Dimension(btnW2,0));
      btnPanel42.setLayout(new BorderLayout(gap1,gap1));
      btn42 = new JButton();
      btnPanel42.add(btn42, BorderLayout.CENTER);

      kPanel40.add(colPanel4);
      kPanel40.add(navPanel4);
      kPanel41.add(cbxPanel4, BorderLayout.CENTER);
      kPanel41.add(btnPanel42, BorderLayout.EAST);
      kPanel4.add(kPanel40, BorderLayout.WEST);
      kPanel4.add(txtPanel4, BorderLayout.CENTER);
      kPanel4.add(kPanel41, BorderLayout.EAST);

//......................................
      JPanel kPanel5 = new JPanel(); // 1 row of the key
      JPanel kPanel50 = new JPanel(); // holds colour panel & arrow button
      JPanel kPanel51 = new JPanel(); // holds chkbx and remove button
      kPanel5.setLayout(new BorderLayout(gap,gap));
      kPanel50.setLayout(new BoxLayout(kPanel50, BoxLayout.X_AXIS));
      kPanel51.setLayout(new BorderLayout(gap1,gap1));
      kPanel50.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      kPanel51.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));

      colBtn5 = new JButton();
      colBtn5.setBackground(new Color(_cols[5], hasAlpha));
      JPanel colPanel5 = new JPanel(new GridLayout());
      colPanel5.add(colBtn5);
      colPanel5.setPreferredSize(new Dimension(colW,0));
      JPanel txtPanel5 = new JPanel();
      txtPanel5.setLayout(new BorderLayout(gap1,gap1));
      txtPanel5.setBorder(BorderFactory.createEtchedBorder(
                                          EtchedBorder.LOWERED));
      tf5 = new ScrollableTextField("", 50);
      tf5.setEditable(false);
      tf5.setBackground(notvis);
      txtPanel5.add(tf5, BorderLayout.CENTER);
      JPanel navPanel5 = new JPanel();
      navPanel5.setPreferredSize(new Dimension(btnW3,0));
      navPanel5.setLayout(new BoxLayout(navPanel5, BoxLayout.X_AXIS));
      JPanel btnPanel50 = new JPanel();
      btnPanel50.setPreferredSize(new Dimension(btnW1,0));
      btnPanel50.setLayout(new BorderLayout(gap1,gap1));
      btn50 = new JButton();
      if(leftIcon != null) {
	 btn50.setIcon(leftIcon);
      }
      btnPanel50.add(btn50, BorderLayout.CENTER);
      JPanel btnPanel51 = new JPanel();
      btnPanel51.setPreferredSize(new Dimension(btnW1,0));
      btnPanel51.setLayout(new BorderLayout(gap1,gap1));
      btn51 = new JButton();
      if(leftIcon != null) {
	 btn51.setIcon(rightIcon);
      }
      btnPanel51.add(btn51, BorderLayout.CENTER);
      navPanel5.add(btnPanel50);
      navPanel5.add(btnPanel51);
      JPanel cbxPanel5 = new JPanel();
      cbxPanel5.setPreferredSize(new Dimension(btnW1,0));
      cbxPanel5.setLayout(new FlowLayout(FlowLayout.CENTER, gap1,gap1));
      cbx5 = new JCheckBox();
      cbxPanel5.add(cbx5, BorderLayout.CENTER);
      cbx5.setSelected(cbx5State);
      JPanel btnPanel52 = new JPanel();
      btnPanel52.setPreferredSize(new Dimension(btnW2,0));
      btnPanel52.setLayout(new BorderLayout(gap1,gap1));
      btn52 = new JButton();
      btnPanel52.add(btn52, BorderLayout.CENTER);

      kPanel50.add(colPanel5);
      kPanel50.add(navPanel5);
      kPanel51.add(cbxPanel5, BorderLayout.CENTER);
      kPanel51.add(btnPanel52, BorderLayout.EAST);
      kPanel5.add(kPanel50, BorderLayout.WEST);
      kPanel5.add(txtPanel5, BorderLayout.CENTER);
      kPanel5.add(kPanel51, BorderLayout.EAST);

//......................................
      setToolTips();
//......................................

      kTopPanel.add(kPanel0);
      kTopPanel.add(kPanel1);
      kTopPanel.add(kPanel2);
      kTopPanel.add(kPanel3);
      kTopPanel.add(kPanel4);
      kTopPanel.add(kPanel5);
      topPanel.add(kTopPanel, BorderLayout.CENTER);
   } // makeAnatKeyGUI()

//......................................

   private void setToolTips() {

      colBtn0.setToolTipText(colTipText);
      colBtn1.setToolTipText(colTipText);
      colBtn2.setToolTipText(colTipText);
      colBtn3.setToolTipText(colTipText);
      colBtn4.setToolTipText(colTipText);
      colBtn5.setToolTipText(colTipText);

      btn00.setToolTipText(leftTipText);
      btn10.setToolTipText(leftTipText);
      btn20.setToolTipText(leftTipText);
      btn30.setToolTipText(leftTipText);
      btn40.setToolTipText(leftTipText);
      btn50.setToolTipText(leftTipText);

      btn01.setToolTipText(rightTipText);
      btn11.setToolTipText(rightTipText);
      btn21.setToolTipText(rightTipText);
      btn31.setToolTipText(rightTipText);
      btn41.setToolTipText(rightTipText);
      btn51.setToolTipText(rightTipText);

      cbx0.setToolTipText(visTipText);
      cbx1.setToolTipText(visTipText);
      cbx2.setToolTipText(visTipText);
      cbx3.setToolTipText(visTipText);
      cbx4.setToolTipText(visTipText);
      cbx5.setToolTipText(visTipText);

      btn02.setToolTipText(zapTipText);
      btn12.setToolTipText(zapTipText);
      btn22.setToolTipText(zapTipText);
      btn32.setToolTipText(zapTipText);
      btn42.setToolTipText(zapTipText);
      btn52.setToolTipText(zapTipText);

   }

//......................................

   public void setZapToolTips(int indx, boolean zapped) {

      if(zapped) {
         zapTipText = "replace in key";
      } else {
         zapTipText = "remove from key";
      }

      switch(indx) {
      case 0:
	 btn02.setToolTipText(zapTipText);
	 break;
      case 1:
	 btn12.setToolTipText(zapTipText);
	 break;
      case 2:
	 btn22.setToolTipText(zapTipText);
	 break;
      case 3:
	 btn32.setToolTipText(zapTipText);
	 break;
      case 4:
	 btn42.setToolTipText(zapTipText);
	 break;
      case 5:
	 btn52.setToolTipText(zapTipText);
	 break;
      default:
         break;
      }

   }
//-------------------------------------------------------------
// inner classes for event handling
//-------------------------------------------------------------
   protected abstract class ButtonHandler1 implements ActionListener {
   }

   protected abstract class ButtonHandler2 implements ActionListener {
   }

   protected abstract class ChBoxHandler implements ItemListener {
   }

} // class AnatKeyGUI
