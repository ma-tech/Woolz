package sectionViewer;

import java.util.*;

import java.awt.*;
import java.awt.event.ComponentListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ActionEvent;
import javax.swing.*;

/**
 *   A table which indicates the colour of anatomy components in SectionViewers.
 *   AnatKey has controls to allow the colour to be changed
 *   and the component visibility to be toggled off / on.
 *   AnatKey is a singleton class.
 */
public class AnatKey extends AnatKeyGUI{

   /**
    *   The one and only instance of AnatKey (per parent).
    */
   protected AnatKey _instance = null;

   /**
    *   The collection of rows in the AnatKey.
    */
   protected  Vector<KeyEntry> _keyEntryVec = null;

   /**
    *   The number of rows in the AnatKey.
    */
   protected  int _nrows = 0;

   /**
    *   A unique index number for KeyElements.
    */
   protected  int _indx = 0;

   protected int _maxW = 600;
   protected int _maxH = 371;
   protected int _minW = 100;
   protected int _minH = 20;

   protected int _userW = _goldenW;
   protected int _userH = _goldenH;
   protected Point _userLoc = new Point(0,0);

//-------------------------------------------------------------
   /**
    *   Constructs a 2D AnatKey with the default title.
    *   Only 1 instance is allowed per parent.
    *   It is the responsibility of the parent to ensure
    *   that only 1 instance is constructed. !!!
    */
   public AnatKey() {
      this(false);
   }

//............................
   /**
    *   Constructs an AnatKey with the default title
    *   and 2D / 3D as specified.
    *   It is the responsibility of the parent to ensure
    *   that only 1 instance is constructed. !!!
    */
   public AnatKey(boolean is3D) {

      super("Anatomy Key", is3D);
      _is3D = is3D;
      _keyEntryVec = new Vector<KeyEntry>();
      ResizeListener resizeListen = new ResizeListener(this);
      this.addComponentListener(resizeListen);
   }

//-------------------------------------------------------------
   /**
    *   Returns the data model for AnatKey.
    *   @return _keyEntryVec
    */
    public Vector<KeyEntry> getModel() {
       return _keyEntryVec;
    }

//-------------------------------------------------------------
   /**
    *   Sets the data model for AnatKey.
    *   @param entries
    */
    public void setModel(Vector<KeyEntry> entries) {
       _keyEntryVec = entries;
    }

//-------------------------------------------------------------
   /**
    *
    */
   public void addRow(String txt) {

      //System.out.println("enter addRow++++++++++++++++++++++++++++++++++++++");

      clear();

      KeyEntry row = new KeyEntry(_indx, _is3D);
      row.setText(txt);
      row.setCol(nextCol());
      _keyEntryVec.add(row);
      _nrows = _keyEntryVec.size();

      if(_nrows == 1) {
         kScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
      }
      
      for(int i=0; i<_nrows; i++) {
         row = (KeyEntry)_keyEntryVec.elementAt(i);
	 kTopPanel.add(row);
      }

      this.pack();
      setUserSize(_userW, _userH);
      setUserLocation(_userLoc);

      //System.out.println("exit addRow");
   }
//-------------------------------------------------------------
   /**
    *
    */
   public void removeRow(int indx) {

      //System.out.println("enter removeRow -------------------------------");

      /* find row with appropriate indx from _keyEntryVec */
      KeyEntry row = getRow(indx);
      _keyEntryVec.remove(row);
      kTopPanel.remove(row);

      _nrows--;
      if(_nrows <= 0) {
         kScroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_NEVER);
      }

      this.pack();
      setUserSize(_userW, _userH);
      setUserLocation(_userLoc);

      //System.out.println("exit removeRow");
   }
//-------------------------------------------------------------
   /**
    *   Returns the number of rows in the anatomy key.
    *   @return The int number of rows in the key.
    */
   public int getNRows() {
      return _nrows;
   }

//-------------------------------------------------------------
   /**
    *   Returns the anatomy key header.
    *   @return The header.
    */
   public void addHeader() {
      if(_header == null) {
         _header = new KeyHeader();
      }
      topPanel.add(_header);
   }

//-------------------------------------------------------------
   /**
    *   Removes all entries and the header.
    */
   public void clear() {
      kTopPanel.removeAll();
   }

//-------------------------------------------------------------
   /**
    *   Sets the visibility of a KeyEntry to the given state.
    *   @param indx the index of the KeyEntry to be made visible.
    */
   public void setEntryText(int indx) {

      KeyEntry row = this.getRow(indx);
      row.setEntryText();

   }

//-------------------------------------------------------------
   /**
    *   Sets the initial colour of the anatomy component
    *   to 1 of 6 distinct colours.
    *   @return the next colour.
    */
   public Color nextCol() {

      int i = (_ncols + _indx -1) % _ncols;
      return new Color(_rgbt[i][0],
                       _rgbt[i][1],
                       _rgbt[i][2]);
   }

//-------------------------------------------------------------
   /**
    *   Sets the colour of the anatomy component
    *   to that returned by a standard color chooser dialog.
    *   The color chooser dialog is opened by clicking on the
    *   Colour chooser button in the AnatKey.
    *   @param col the new colour of the anatomy component.
    *   @param indx the position (i.e. the row number in the AnatKey)
    *   of the anatomy component whose colour is being changed.
    */
   public void setCol(int indx) {

      Color col = null;
      KeyEntry row = null;

      row = getRow(indx);

      if(row != null) {
	 col = JColorChooser.showDialog(null,
		     "choose colour for anatomy component",
		     row.getCol());

	 if(col != null) {
	    row.setCol(col);
	 }
      }

   }

//-------------------------------------------------------------
   /**
    *   Returns the colour of an anatomy component at
    *   the given index in the Anatomy Key.
    *   @param indx the index.
    *   @return a new Color object.
    */
   public Color getColor(int indx){

      Color col = null;
      KeyEntry row = null;

      row = getRow(indx);
      if(row != null) {
	col = row.getCol();
      }

      return col;
   }

//-------------------------------------------------------------
   /**
    *   Resets the AnatKey to its empty state.
    */
   public void reset() {

      int num = _keyEntryVec.size();
      int indxs[] = new int[num];
      KeyEntry row = null;

      /* first get the collection of index numbers */
      for(int i = 0; i<num; i++) {
         row = (KeyEntry)_keyEntryVec.elementAt(i);
	 indxs[i] = row.getIndx();
      }
      row = null;

      /* then remove each row */
      for(int j = 0; j<num; j++) {
         removeRow(indxs[j]);
      }
   }

//-------------------------------------------------------------
   /**
    *   Returns the keyEntry with the specified indx.
    *   @param indx a unique index for each keyEntry.
    *   @return the KeyEntry with the specified indx.
    */
   public KeyEntry getRow(int indx) {

      int i = -1;
      int size = _keyEntryVec.size();
      for(i=0; i<size; i++) {
	 if(((KeyEntry)_keyEntryVec.elementAt(i)).getIndx() == indx) {
	    break;
	 }
      }
      return (KeyEntry)_keyEntryVec.elementAt(i);
   }

//-------------------------------------------------------------
   /**
    *   Increments the unique index and returns the new value.
    *   @return _indx.
    */
   public int nextIndx() {
      _indx++;
      return _indx;
   }

//-------------------------------------------------------------
   /**
    *   Sets all the 3D visibility control icons.
    *   @param state true if controls are to have viz3DIcon.
    */
   public void set3DVisIcons(boolean visible) {
      int size = _keyEntryVec.size();
      for(int i=0; i<size; i++) {
	 ((KeyEntry)_keyEntryVec.elementAt(i)).set3DVisIcon(visible);
      }
   }

//-------------------------------------------------------------
   private void setUserLocation(Point p) {
      this.setLocation(p.x, p.y);
   }

   private void setUserSize(int w, int h) {
      this.setPreferredSize(new Dimension(w, h));
   }
//-------------------------------------------------------------
   /**
    *   For testing only.
    */
/*
  public void main(String argv[]) {

	AnatKey _key;

        _key = AnatKey.instance();
	_key.setSize(300,500);
	_key.pack();
	_key.setVisible(true);

  }
*/
   /**
    *   Listens for re-size and move events, but only every once in a while.
    *   The timer is used to set the interval between events that are actioned.
    *   Resizing should only throw an event when the resize has finished,
    *   but on Linux it throws events continually.
    *   Moving throws events continually on all platforms.:w
    */
   class ResizeListener implements ComponentListener {

      private javax.swing.Timer timer;
      private AnatKey key;
      private static final int DELAY = 200;

      public ResizeListener(AnatKey key) {
         this.key = key;

	 timer = new javax.swing.Timer(DELAY, new AbstractAction() {
	    public void actionPerformed(ActionEvent e) {
	       savePrefs();
	    }
	 });
	 timer.setRepeats(false);
      }

      /**
       * @param e
       */
      public void componentResized(ComponentEvent e) {
	 timer.start();
      }

      private void savePrefs() {
         if(!key.isShowing()) {
	    return;
	 }
	 _userW = key.getWidth();
	 _userH = key.getHeight();
	 _userLoc = key.getLocationOnScreen();
      }

      public void componentMoved(ComponentEvent e) {
	 //System.out.println("move");
	 timer.start();
      }
      public void componentShown(ComponentEvent e) {
      }
      public void componentHidden(ComponentEvent e) {
      }
   } // class ResizeListener

} // class AnatKey
