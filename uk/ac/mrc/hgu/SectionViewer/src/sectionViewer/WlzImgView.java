package sectionViewer;
import sectionViewer.*;

import java.awt.*;
import java.awt.color.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.io.*;
import java.util.*;

import wsetter.*;
import uk.ac.mrc.hgu.Wlz.*;

//-------------------------------------------------------------
// 'view' class for wlz _objects
// uses RGB colour model (instead of indexed)
//-------------------------------------------------------------

public class WlzImgView extends Component {

   private final boolean _debug = false;

   private WlzObject _obj = null;
   private WlzIBox2 _bBox = null;
   private WlzIVertex2 _org = null;
   private WlzIVertex2 _oorg = null;
   private WlzIVertex2 _torg = null;
   private WlzIVertex2 _aorg[] = null;
   private WlzGreyValueWSpace _gVWSp = null;
   private BufferedImage _bufImage = null;
   private BufferedImage _obufImage = null;
   private BufferedImage _tbufImage = null;
   private BufferedImage _abufImage[] = null;

   private double _mag;
   private double _ixofs;
   private double _iyofs;
   private double _oxofs;
   private double _oyofs;
   private double _txofs;
   private double _tyofs;
   private double _axofs[] = null;
   private double _ayofs[] = null;
   private double _xofsGL;
   private double _yofsGL;
   private double _xofsTC;
   private double _yofsTC;
   private double _xofsFP;
   private double _yofsFP;
   private double _xofsAP;
   private double _yofsAP;

   private String _objStats;
   private int _greyVal;
   private Point _pos = null;
   private Rectangle _bRec = null;
   private Rectangle _obRec = null;
   private Rectangle _tbRec = null;
   private Rectangle _abRec[] = null;

   private Dimension _imgSize = null;

   private DirectColorModel _dcm = null;
   private WritableRaster _ras = null;
   private WritableRaster _oras = null;
   private WritableRaster _tras = null;
   private WritableRaster _aras[] = null;
   private DataBuffer _dataBuf = null;
   private DataBuffer _odataBuf = null;
   private DataBuffer _tdataBuf = null;
   private DataBuffer _adataBuf[] = null;
   private SinglePixelPackedSampleModel _sppsm = null;
   private SinglePixelPackedSampleModel _osppsm = null;
   private SinglePixelPackedSampleModel _tsppsm = null;
   private SinglePixelPackedSampleModel _asppsm[] = null;
   private int _nbits = 32;
   private int _masks[] = {0xff0000, 0xff00, 0xff, 0xff000000}; // (R G B alpha)
   private int _imageData[] = null;
   private int _oimageData[] = null;
   private int _timageData[] = null;
   private int _aimageData[][] = null;

   private GeneralPath _threshConstraintPath = null;
   private int _npts = 0;
   private int _num;
   private int _fpr = 4; // fixed point indicator radius
   private int _apr = 4; // axis point indicator radius

   private boolean _viz[]; // visibility of anatomy objects

   private Vector _intersectionVec = null;
   private Vector _interColVec = null;
   private Vector _fixedPointVec = null;
   private Vector _axisPointVec = null;

   private boolean _fixedPoint;
   private boolean _axisPoint;
   private boolean _axis;
   private boolean _intersection;
   private boolean _anatomy;
   private boolean _threshold;
   private boolean _threshConstraint;
   private boolean _overlay;
   private boolean _showInterSecLines = true;

   private boolean _tiepoint = true;
   private Vector tps = null;
   private Vector tpsCol = null;

   //-------------------------------------------------------------
   // constructor
   public WlzImgView() {

      _fixedPoint = false;
      _axisPoint = false;
      _axis = false;
      _intersection = false;
      _anatomy = false;
      _threshold = false;
      _threshConstraint = false;
      _overlay = false;

      _objStats = "";
      _greyVal = -1;
      _pos = new Point();
      _mag = 1.0;
      _imgSize = new Dimension();

      _num = AnatKey._nrows;

      _abufImage = new BufferedImage[_num];
      _adataBuf = new DataBuffer[_num];
      _asppsm = new SinglePixelPackedSampleModel[_num];
      _aimageData = new int[_num][];
      _aorg = new WlzIVertex2[_num];
      _axofs = new double[_num];
      _ayofs = new double[_num];
      _abRec = new Rectangle[_num];
      _aras = new WritableRaster[_num];

      // initialise
      for(int i=0; i<_num; i++) {
	 _abufImage[i] = null;
	 _adataBuf[i] = null;
	 _asppsm[i] = null;
	 _aimageData[i] = null;
	 _aorg[i] = null;
	 _axofs[i] = 0.0;
	 _ayofs[i] = 0.0;
	 _abRec[i] = null;
	 _aras[i] = null;
      }
      _intersectionVec = new Vector();
      _interColVec = new Vector();
      _interColVec.add(Color.red);
      _fixedPointVec = new Vector();
      _axisPointVec = new Vector();

   }

   //-------------------------------------------------------------
   // makes buffered image from the 2D section
   public void setWlzObj(WlzObject newObj) throws WlzException {
      _obj = newObj;
      _bBox = WlzObject.WlzBoundingBox2I(_obj);

      if(_debug == true) System.out.println("setWlzObj");
      /*
	 System.out.println("xmin = "+_bBox.xMin);
	 System.out.println("xmax = "+_bBox.xMax);
	 System.out.println("ymin = "+_bBox.yMin);
	 System.out.println("ymax = "+_bBox.yMax);
       */

      //_org = new WlzIVertex2(0,0);
      _org = new WlzIVertex2(_bBox.xMin, _bBox.yMin);
      _gVWSp = WlzObject.WlzGreyValueMakeWSp(newObj);
      _bRec = _bBox.toRectangle();
      /*
	 System.out.println("w = "+_bRec.width);
	 System.out.println("h = "+_bRec.height);
       */
      _imgSize.setSize(_bRec.width, _bRec.height);
      setSize(_imgSize);


      _dcm = new DirectColorModel(
	    _nbits,
	    _masks[0],
	    _masks[1],
	    _masks[2],
	    _masks[3]);
      //System.out.println(_dcm.toString());

      _imageData = new int[_bRec.width * _bRec.height];
      byte dstArrayDat[][][] = new byte[1][][] ;
      dstArrayDat[0] = null;
      WlzIVertex2 size = new WlzIVertex2(_bBox.xMax - _bBox.xMin + 1,
	    _bBox.yMax - _bBox.yMin + 1);
      WlzIVertex2 dstSize[] = new WlzIVertex2[1];
      dstSize[0] = null;
      WlzObject.WlzToUArray2D(dstSize, dstArrayDat, _obj, _org, size, 0);

      // remember these are grey level images
      byte bb;
      int ii;
      for(int idY = 0; idY < _bRec.height; ++idY) {
	 for(int idX = 0; idX < _bRec.width; ++idX) {

	    // convert byte (8 bits) to 2's complement int (32 bits)
	    bb = (byte)(dstArrayDat[0][idY][idX]);
	    ii = bb & 0xff;
	    // copy same value to red green and blue components
	    // and set alpha (255 => opaque)
	    _imageData[(idY * _bRec.width) + idX] =
	       ii | (ii << 8) | (ii << 16) | (255 << 24);
	 } // for
      } // for

      makeBufImage();
   }

   //-------------------------------------------------------------
   public void makeBufImage() {

      _dataBuf = new DataBufferInt(_imageData,
	    _bRec.width * _bRec.height);
      if(_dataBuf == null) {
	 System.out.println("_dataBuf is null");
      }

      _sppsm = new SinglePixelPackedSampleModel(
	    DataBuffer.TYPE_INT,
	    _bRec.width,
	    _bRec.height,
	    _masks);

      _ras = Raster.createWritableRaster(_sppsm,
	    _dataBuf,
	    null);

      try {
	 _bufImage = new BufferedImage(_dcm, _ras, false, null);
      }
      catch(IllegalArgumentException e) {
	 System.out.println(e.getMessage());
      }
      if(_bufImage == null) {
	 System.out.println("_bufImage is null");
      }

      setGreyLevelOffsets();
   }

   //-------------------------------------------------------------
   // makes image raster from overlay 2D section
   // which should have the same view structure as wlz _object
   public void setOverlayObj(WlzObject newObj) throws WlzException {
      WlzObject obj = newObj;
      WlzIBox2 bBox = null;
      Dimension imgSize = new Dimension();
      int imageData[] = null;
      byte dstArrayDat[][][] = new byte[1][][] ;
      dstArrayDat[0] = null;

      try {
	 bBox = WlzObject.WlzBoundingBox2I(obj);
	 _oorg = new WlzIVertex2(bBox.xMin, bBox.yMin);
	 _obRec = bBox.toRectangle();

	 _oimageData = new int[_obRec.width * _obRec.height];
	 WlzIVertex2 size = new WlzIVertex2(bBox.xMax - bBox.xMin + 1,
	       bBox.yMax - bBox.yMin + 1);
	 WlzIVertex2 dstSize[] = new WlzIVertex2[1];
	 dstSize[0] = null;
	 WlzObject.WlzToUArray2D(dstSize, dstArrayDat, obj, _oorg, size, 0);

	 byte bb;
	 int ii;
	 for(int idY = 0; idY < _obRec.height; ++idY) {
	    for(int idX = 0; idX < _obRec.width; ++idX) {

	       // convert byte (8 bits) to 2's complement int (32 bits)
	       bb = (byte)(dstArrayDat[0][idY][idX]);
	       ii = bb & 0xff;
	       if(ii == 0) { // change black background to white
		  _oimageData[(idY * _obRec.width) + idX] =
		     255 | (255 << 8) | (255 << 16) | (0 << 24);
	       } else {
		  // colour the overlay part
		  _oimageData[(idY * _obRec.width) + idX] =
		     //0 | (255 << 8) | (0 << 16) | (125 << 24); // green
		     //0 | (255 << 8) | (255 << 16) | (125 << 24); // yellow
		     //255 | (255 << 8) | (0 << 16) | (125 << 24); // cyan
		     255 | (0 << 8) | (255 << 16) | (125 << 24); // magenta
	       }
	    } // for
	 } // for

	 makeBufOImage();

      } // try
      catch (WlzException e) {
	 // gives an error if the anatomy doesn't exist
	 //System.out.println(e.getMessage());
      }
   }

   //-------------------------------------------------------------
   public void makeBufOImage() {

      // buffered image for overlay
      _odataBuf = new DataBufferInt(_oimageData,
	    _obRec.width * _obRec.height);
      if(_odataBuf == null) {
	 System.out.println("_odataBuf is null");
      }

      _osppsm = new SinglePixelPackedSampleModel(
	    DataBuffer.TYPE_INT,
	    _obRec.width,
	    _obRec.height,
	    _masks);

      _oras = Raster.createWritableRaster(_osppsm,
	    _odataBuf,
	    null);

      try {
	 _obufImage = new BufferedImage(_dcm, _oras, false, null);
      }
      catch(IllegalArgumentException e) {
	 System.out.println(e.getMessage());
      }
      if(_obufImage == null) {
	 System.out.println("_obufImage is null");
      }

      setOverlayOffsets();
      _overlay = true;
   }

   //-------------------------------------------------------------
   // makes image raster from thresholded 2D section
   // which should have the same view structure as wlz _object
   public void setThresholdObj(WlzObject newObj) throws WlzException {
      WlzObject obj = newObj;
      WlzIBox2 bBox = null;
      Dimension imgSize = new Dimension();
      int imageData[] = null;
      byte dstArrayDat[][][] = new byte[1][][] ;
      dstArrayDat[0] = null;

      try {
	 bBox = WlzObject.WlzBoundingBox2I(obj);
	 /*
         System.out.println("thresholded obj ...");
         System.out.println("xmin: "+bBox.xMin+", xmax: "+bBox.xMax);
         System.out.println("ymin: "+bBox.yMin+", ymax: "+bBox.yMax);
	 */
	 _torg = new WlzIVertex2(bBox.xMin, bBox.yMin);
	 _tbRec = bBox.toRectangle();

	 _timageData = new int[_tbRec.width * _tbRec.height];
	 WlzIVertex2 size = new WlzIVertex2(bBox.xMax - bBox.xMin + 1,
	       bBox.yMax - bBox.yMin + 1);
	 WlzIVertex2 dstSize[] = new WlzIVertex2[1];
	 dstSize[0] = null;
	 WlzObject.WlzToUArray2D(dstSize, dstArrayDat, obj, _torg, size, 0);

	 byte bb;
	 int ii;
	 for(int idY = 0; idY < _tbRec.height; ++idY) {
	    for(int idX = 0; idX < _tbRec.width; ++idX) {

	       // convert byte (8 bits) to 2's complement int (32 bits)
	       bb = (byte)(dstArrayDat[0][idY][idX]);
	       ii = bb & 0xff;
	       if(ii == 255) { // ignore white (background?)
		  _timageData[(idY * _tbRec.width) + idX] =
		     255 | (255 << 8) | (255 << 16) | (0 << 24);
	       } else {
		  // colour the thresholded part
		  _timageData[(idY * _tbRec.width) + idX] =
		     0 | (255 << 8) | (0 << 16) | (125 << 24); // green
		  // | (255 << 8) | (255 << 16) | (125 << 24); // yellow
		  //255 | (255 << 8) | (0 << 16) | (125 << 24); // cyan
		  //255 | (0 << 8) | (255 << 16) | (125 << 24); // magenta
	       }
	    } // for
	 } // for

	 makeBufTImage();

      } // try
      catch (WlzException e) {
	 //System.out.println("Wlz error #A");
	 //System.out.println(e.getMessage());
      }
   }

   //-------------------------------------------------------------
   public void makeBufTImage() {

      // make buffered image for thresholded image
      _tdataBuf = new DataBufferInt(_timageData,
	    _tbRec.width * _tbRec.height);
      if(_tdataBuf == null) {
	 System.out.println("_tdataBuf is null");
      }

      _tsppsm = new SinglePixelPackedSampleModel(
	    DataBuffer.TYPE_INT,
	    _tbRec.width,
	    _tbRec.height,
	    _masks);

      _tras = Raster.createWritableRaster(_tsppsm,
	    _tdataBuf,
	    null);
      try {
	 _tbufImage = new BufferedImage(_dcm, _tras, false, null);
      }
      catch(IllegalArgumentException e) {
	 System.out.println(e.getMessage());
      }
      if(_tbufImage == null) {
	 System.out.println("_tbufImage is null");
      }

      setThresholdOffsets();
      _threshold = true;
   }

   //-------------------------------------------------------------
   // makes image raster(s) from anatomy section(s)
   // which should have the same view structure as wlz _object
   public void setAnatomyObj(WlzObject[] objs, boolean[] viz)
      throws WlzException {

	 int total = AnatKey._nrows;
	 _viz = viz;

	 for(int i = 0; i < total; i++) {
	    if(_viz[i] == false) {
	       _aorg[i] = null;
	       _abufImage[i] = null;
	       continue;
	    }
	    WlzObject obj = objs[i];
	    if(obj != null) {
	       WlzIBox2 bBox = null;
	       Dimension imgSize = new Dimension();
	       int imageData[] = null;
	       byte dstArrayDat[][][] = new byte[1][][] ;
	       dstArrayDat[0] = null;

	       try {
		  bBox = WlzObject.WlzBoundingBox2I(obj);
		  _aorg[i] = new WlzIVertex2(bBox.xMin, bBox.yMin);
		  _abRec[i] = bBox.toRectangle();

		  _aimageData[i] = new int[_abRec[i].width *
		     _abRec[i].height];
		  WlzIVertex2 size = new WlzIVertex2(bBox.xMax - bBox.xMin + 1,
			bBox.yMax - bBox.yMin + 1);
		  WlzIVertex2 dstSize[] = new WlzIVertex2[1];
		  dstSize[0] = null;
		  WlzObject.WlzToUArray2D(dstSize,
			dstArrayDat,
			obj,
			_aorg[i],
			size,
			0);

		  byte bb;
		  int ii;
		  for(int idY = 0; idY < _abRec[i].height; ++idY) {
		     for(int idX = 0; idX < _abRec[i].width; ++idX) {

			// convert byte (8 bits) to 2's complement int (32 bits)
			bb = (byte)(dstArrayDat[0][idY][idX]);
			ii = bb & 0xff;
			if(ii == 0) { // change black background to white
			   _aimageData[i][(idY * _abRec[i].width) + idX] =
			      255 | (255 << 8) | (255 << 16) | (0 << 24);
			} else {
			   // colour the overlay part
			   _aimageData[i][(idY * _abRec[i].width) + idX] =
			      AnatKey._cols[i];
			}
		     } // for
		  } // for

		  makeBufAImage(i);

	       } // try
	       catch (WlzException e) {
		  /*
		  //thrown if anatomy isn't intersected by section
		  System.out.println("setAnatomy ...");
		  System.out.println(e.getMessage());
		   */
	       }
	    } else {
	       _aorg[i] = null;
	       _abufImage[i] = null;
	    } // if(obj != null)

	 } // for

      } // setAnatomyObj()

   //-------------------------------------------------------------
   public void makeBufAImage(int indx) {

      // buffered image for overlay
      _adataBuf[indx] = new DataBufferInt(_aimageData[indx],
	    _abRec[indx].width * _abRec[indx].height);
      if(_adataBuf[indx] == null) {
	 System.out.println("_adataBuf["+indx+"] is null");
      }

      _asppsm[indx] = new SinglePixelPackedSampleModel(
	    DataBuffer.TYPE_INT,
	    _abRec[indx].width,
	    _abRec[indx].height,
	    _masks);

      _aras[indx] = Raster.createWritableRaster(_asppsm[indx],
	    _adataBuf[indx],
	    null);

      try {
	 _abufImage[indx] = new BufferedImage(_dcm, _aras[indx], false, null);
      }
      catch(IllegalArgumentException e) {
	 System.out.println(e.getMessage());
      }
      if(_abufImage[indx] == null) {
	 System.out.println("_abufImage is null");
      }

      setAnatomyOffsets();
      _anatomy = true;

   }

   //-------------------------------------------------------------
   // makes a GeneralPath from Vector of Floats
   public void setThreshConstraint(Vector xPoints, Vector yPoints) {

      int npts = 0;

      npts = (yPoints.size() >= xPoints.size()) ?
                      xPoints.size() : yPoints.size();


      if(_threshConstraintPath == null) {
	 _threshConstraintPath = new GeneralPath(
	                  GeneralPath.WIND_NON_ZERO, npts);
         _threshConstraintPath.moveTo(
	     ((Float)xPoints.elementAt(0)).floatValue(),
	     ((Float)yPoints.elementAt(0)).floatValue());
	 _npts = 1;

      }


      for(int i=_npts; i<npts; i++) {
         _threshConstraintPath.lineTo(
	     ((Float)xPoints.elementAt(i)).floatValue(),
	     ((Float)yPoints.elementAt(i)).floatValue());
      }

      _npts = npts;

   }
   //-------------------------------------------------------------
   public void closeThreshConstraint() {
      _threshConstraintPath.closePath();
   }
   //-------------------------------------------------------------
   public void getPolyDomain(GeneralPath GP) {

   }
   //-------------------------------------------------------------
   public void paint(Graphics g) {
      Graphics2D g2 = (Graphics2D)g;
      g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER));
      g2.scale(_mag, _mag);
      drawGreyImage(g2);
      drawOverlay(g2);
      drawAnatomy(g2);
      drawIntersection(g2);
      drawFixedPoint(g2);
      drawAxisPoint(g2);
      drawAxis(g2);
      // temporarily disabled ... don't remove
      drawThreshold(g2);
      drawThreshConstraint(g2);
      //--------------------------------------
      drawTiePoint(g2);
   }

   //-------------------------------------------------------------
   public void drawGreyImage(Graphics2D g2) {
      if(_bufImage != null) {
	 g2.translate(_xofsGL, _yofsGL);
	 g2.drawImage(_bufImage, 0, 0, null);
	 g2.translate(-_xofsGL, -_yofsGL);
      }
   }

   //-------------------------------------------------------------
   public void drawOverlay(Graphics2D g2) {
      if (_overlay == false) return;
      if(_obufImage != null) {
	 g2.translate(_oxofs, _oyofs);
	 g2.drawImage(_obufImage, 0, 0, null);
	 g2.translate(-_oxofs, -_oyofs);
      }
   }

   //-------------------------------------------------------------
   public void drawThreshold(Graphics2D g2) {
      if (_threshold == false) return;
      if(_tbufImage != null) {
	 g2.translate(_txofs, _tyofs);
	 g2.drawImage(_tbufImage, 0, 0, null);
	 g2.translate(-_txofs, -_tyofs);
      }
   }

   //-------------------------------------------------------------
   public void drawAnatomy(Graphics2D g2) {
      if (_anatomy == false) return;
      int total = AnatKey._nrows;

      for(int i=0; i<total; i++) {
	 if(_viz[i] == false) continue;
	 if(_abufImage[i] != null) {
	    g2.translate(_axofs[i], _ayofs[i]);
	    g2.drawImage(_abufImage[i], 0, 0, null);
	    g2.translate(-_axofs[i], -_ayofs[i]);
	 }
      }
   }

   //-------------------------------------------------------------
   public void drawIntersection(Graphics2D g2) {

      if(_intersection == false) return;
      int num = _intersectionVec.size();
      if(num == 0) return;
      boolean COL = false;
      int numCols = _interColVec.size();
      if(numCols >= num) COL = true;

      Line2D.Double line = null;

      g2.translate(_ixofs, _iyofs);
      for(int i=0; i<num; i++) {
	 line = (Line2D.Double)_intersectionVec.elementAt(i);
	 if(line == null) { // maybe two views are parallel
	    g2.translate(-_ixofs, -_iyofs);
	    return;
	 }
	 if(COL) g2.setColor((Color)_interColVec.elementAt(i));
	 g2.drawLine((int)line.getX1(),
	       (int)line.getY1(),
	       (int)line.getX2(),
	       (int)line.getY2());
      }
      g2.translate(-_ixofs, -_iyofs);
      if(COL) g2.setColor(Color.black);
   }

   //-------------------------------------------------------------
   public void drawFixedPoint(Graphics2D g2) {

      if(_fixedPoint == false) return;
      // draw a cross centred on fixed point
      int num = _fixedPointVec.size();
      //System.out.println("_fixedPoint vec: size = "+num);
      if(num == 0) return;

      Color colOrig = g2.getColor();
      Line2D.Double line = null;

      g2.setColor(Color.green);
      g2.translate(_xofsFP, _yofsFP);
      for(int i=0; i<num; i++) {
	 line = (Line2D.Double)_fixedPointVec.elementAt(i);
	 g2.drawLine((int)line.getX1(),
	       (int)line.getY1(),
	       (int)line.getX2(),
	       (int)line.getY2());
	 if(i == 0) {
	    int x = (int)line.getX1() - _fpr;
	    int y = (int)line.getY1() - _fpr;
	    g2.drawOval(x,y,2*_fpr,2*_fpr);
	 }
      }
      g2.translate(-_xofsFP, -_yofsFP);
      g2.setColor(colOrig);
   }

   //-------------------------------------------------------------
   public void drawAxisPoint(Graphics2D g2) {

      // draw a cross centred on axis point
      if(_axisPoint == false) return;
      int num = _axisPointVec.size();
      //System.out.println("_axisPoint vec: size = "+num);
      if(num == 0) return;

      Color colOrig = g2.getColor();
      Line2D.Double line = null;

      g2.setColor(Color.cyan);
      g2.translate(_xofsAP, _yofsAP);
      for(int i=0; i<num; i++) {
	 line = (Line2D.Double)_axisPointVec.elementAt(i);
	 g2.drawLine((int)line.getX1(),
	       (int)line.getY1(),
	       (int)line.getX2(),
	       (int)line.getY2());
	 if(i == 0) {
	    int x = (int)line.getX1() - _apr;
	    int y = (int)line.getY1() - _apr;
	    g2.drawOval(x,y,2*_apr,2*_apr);
	 }
      }
      g2.translate(-_xofsAP, -_yofsAP);
      g2.setColor(colOrig);
   }

   //-------------------------------------------------------------
   public void drawAxis(Graphics2D g2) {

      // draw a line between fixed point and axis point
      if(_axis == false) return;
      if((_fixedPointVec == null) || (_fixedPointVec.size() == 0)) return;
      if((_axisPointVec == null) || (_axisPointVec.size() == 0)) return;

      Color colOrig = g2.getColor();
      Line2D.Double line = null;

      g2.setColor(Color.cyan);
      g2.translate(_xofsAP, _yofsAP);

      line = (Line2D.Double)_axisPointVec.elementAt(0);
      int x1 = (int)line.getX1();
      int y1 = (int)line.getY1();

      line = (Line2D.Double)_axisPointVec.elementAt(1);

      line = (Line2D.Double)_fixedPointVec.elementAt(0);
      int x2 = (int)line.getX1();
      int y2 = (int)line.getY1();

      line = (Line2D.Double)_fixedPointVec.elementAt(1);

      g2.drawLine(x1,y1,x2,y2);

      g2.translate(-_xofsAP, -_yofsAP);
      g2.setColor(colOrig);
   }

   //-------------------------------------------------------------
   public void drawThreshConstraint(Graphics2D g2) {
      if(_threshConstraint && (_threshConstraintPath != null)) {
	 g2.setColor(Color.red);
	 g2.draw(_threshConstraintPath);
      }
   }

   //-------------------------------------------------------------
   public void clearOverlay() {
      _obufImage = null;
      _overlay = false;
   }

   //-------------------------------------------------------------
   public void enableFixedPoint(boolean state) {
      _fixedPoint = state;
   }

   //-------------------------------------------------------------
   public void enableAxisPoint(boolean state) {
      _axisPoint = state;
   }

   //-------------------------------------------------------------
   public void enableAxis(boolean state) {
      _axis = state;
   }

   //-------------------------------------------------------------
   public void clearIntersection() {
      _intersection = false;
   }

   //-------------------------------------------------------------
   public void clearThreshold() {
      _tbufImage = null;
      _threshold = false;
   }

   //-------------------------------------------------------------
   public void enableThreshConstraint(boolean state) {
      _threshConstraint = state;
   }

   //-------------------------------------------------------------
   public void clearThreshConstraint() {
      _threshConstraintPath = null;
   }

   //-------------------------------------------------------------
   public void clearAnatomy() {
      for(int i=0; i<_num; i++) {
	 _abufImage[i] = null;
      }
      _anatomy = false;
   }

   //-------------------------------------------------------------
   public GeneralPath getThreshConstraint() {
      return _threshConstraintPath;
   }

   //-------------------------------------------------------------
   public void updateStats(double relX, double relY) {

      if((_obj == null) || (_gVWSp == null)) {
	 _objStats = "";
      } else {
	 _pos.setLocation(_bBox.xMin + relX, _bBox.yMin + relY);
	 _greyVal = WlzObject.WlzGreyValueGetI(_gVWSp,
	       0.0, (double )_pos.y, (double )_pos.x);

	 _objStats = "G(" + _pos.x + ", " + _pos.y + ") = " +_greyVal;
      }
      //System.out.println("stats: "+_objStats);

   } // updateStats

   //-------------------------------------------------------------
   public Point getOffSet(){
     return new Point(_bBox.xMin, _bBox.yMin);
   }

   //-------------------------------------------------------------
   protected void setIntersectionVec(Line2D.Double[] lines) {

      if(_intersectionVec == null) {
	 _intersectionVec = new Vector();
      } else {
	 _intersectionVec.clear();
      }

      if((lines == null) || (lines.length == 0)) return;

      for(int i=0; i<lines.length; i++) {
	 _intersectionVec.add(lines[i]);
      }
      setIntersectionOffsets();
      _intersection = true;
   }

   //-------------------------------------------------------------
   protected void setInterColVec(Color[] cols) {

      if((cols == null) || (cols.length == 0)) return;

      if(_interColVec == null) {
	 _interColVec = new Vector();
      } else {
	 _interColVec.clear();
      }
      for(int i=0; i<cols.length; i++) {
	 _interColVec.add(cols[i]);
      }
   }

   //-------------------------------------------------------------
   void drawTiePoint(Graphics2D g){
     if (!_tiepoint) return;
     if (null == tps) return;

     Color orgColor = g.getColor();
     g.translate(_xofsFP, _yofsFP);
     g.scale(1/_mag, 1/_mag);
     for (int i = 0; i < tps.size(); i++){
       Point p = (Point) tps.get(i);
       g.setColor((Color)tpsCol.get(i));
       g.fillOval((int)_mag*p.x - 4, (int)_mag*p.y - 4, 8, 8);
     }
     g.scale(_mag, _mag);
     g.setColor(orgColor);
     g.translate(-_xofsFP, -_yofsFP);
   }

   //-------------------------------------------------------------
   public void setTiePoint(Vector tps, Vector tpsCol){
     this.tps = tps;
     this.tpsCol = tpsCol;
     setFixedPointOffsets();
   }

   //-------------------------------------------------------------
   public void enableTiePoint(boolean state) {
         _tiepoint = state;
   }
   //-------------------------------------------------------------
   protected void setFixedPointVec(double[] fpa) {

      if(_debug) System.out.println("setFixedPointVec");
      Line2D.Double line = null;

      if((fpa == null) || (fpa.length != 3)) return;

      if(_fixedPointVec == null) {
	 _fixedPointVec = new Vector();
      } else {
	 _fixedPointVec.clear();
      }

      line = new Line2D.Double(fpa[0], fpa[1], fpa[0]+_fpr, fpa[1]);
      _fixedPointVec.add(line);
      if(_debug == true) {
	 System.out.print(String.valueOf(line.getX1())+", ");
	 System.out.println(String.valueOf(line.getY1())+", ");
      }
      line = new Line2D.Double(fpa[0], fpa[1], fpa[0], fpa[1]-_fpr);
      _fixedPointVec.add(line);
      line = new Line2D.Double(fpa[0], fpa[1], fpa[0]-_fpr, fpa[1]);
      _fixedPointVec.add(line);
      line = new Line2D.Double(fpa[0], fpa[1], fpa[0], fpa[1]+_fpr);
      _fixedPointVec.add(line);

      setFixedPointOffsets();
      _fixedPoint = true;
   }

   //-------------------------------------------------------------
   protected void setAxisPointVec(double[] apa) {

      if(_debug) System.out.println("setAxisPointVec");
      Line2D.Double line = null;

      if((apa == null) || (apa.length != 3)) return;

      if(_axisPointVec == null) {
	 _axisPointVec = new Vector();
      } else {
	 _axisPointVec.clear();
      }

      line = new Line2D.Double(apa[0], apa[1], apa[0]+_apr, apa[1]);
      _axisPointVec.add(line);
      if(_debug == true) {
	 System.out.print(String.valueOf(line.getX1())+", ");
	 System.out.println(String.valueOf(line.getY1())+", ");
      }
      line = new Line2D.Double(apa[0], apa[1], apa[0], apa[1]-_apr);
      _axisPointVec.add(line);
      line = new Line2D.Double(apa[0], apa[1], apa[0]-_apr, apa[1]);
      _axisPointVec.add(line);
      line = new Line2D.Double(apa[0], apa[1], apa[0], apa[1]+_apr);
      _axisPointVec.add(line);

      setAxisPointOffsets();
      _axisPoint = true;
      _axis = true;
   }

   //-------------------------------------------------------------
   public void printLine(Line2D.Double line) {
      System.out.println(Double.toString(line.getX1())+", "+
                         Double.toString(line.getY1())+", "+
			 Double.toString(line.getX2())+", "+
			 Double.toString(line.getY2()));
   }

   //-------------------------------------------------------------
   public void setGreyLevelOffsets() {

      _xofsGL = 0;
      _yofsGL = 0;

   }
   //-------------------------------------------------------------
   protected void setIntersectionOffsets() {
      _ixofs = 0;
      _iyofs = 0;
   }
   //-------------------------------------------------------------
   public void setFixedPointOffsets() {
      double X = Math.abs((_bBox.xMax - _bBox.xMin));
      double Y = Math.abs((_bBox.yMax - _bBox.yMin));

      _xofsFP = X - _bBox.xMax;
      _yofsFP = Y - _bBox.yMax;
   }
   //-------------------------------------------------------------
   protected double[] getFixedPointOffsets() {
      double ret[] = new double[2];

      ret[0] = _xofsFP;
      ret[1] = _yofsFP;

      return ret;
   }
   //-------------------------------------------------------------
   public void setAxisPointOffsets() {
      double X = Math.abs((_bBox.xMax - _bBox.xMin));
      double Y = Math.abs((_bBox.yMax - _bBox.yMin));

      _xofsAP = X - _bBox.xMax;
      _yofsAP = Y - _bBox.yMax;
   }
   //-------------------------------------------------------------
   protected void setOverlayOffsets() {
      _oxofs = Math.abs(_oorg.vtX - _org.vtX);
      _oyofs = Math.abs(_oorg.vtY - _org.vtY);
   }
   //-------------------------------------------------------------
   protected void setThresholdOffsets() {
      _txofs = Math.abs(_torg.vtX - _org.vtX);
      _tyofs = Math.abs(_torg.vtY - _org.vtY);
   }
   //-------------------------------------------------------------
   protected void setAnatomyOffsets() {
      int total = AnatKey._nrows ;
      if(total==0) return;

      for(int i=0; i<total; i++) {
	 if(_aorg[i] != null) {
	    _axofs[i] = Math.abs(_aorg[i].vtX - _org.vtX);
	    _ayofs[i] = Math.abs(_aorg[i].vtY - _org.vtY);
	 }
      }
   }
   //-------------------------------------------------------------
   public void setThreshConstraintOffsets(double x, double y) {
      _xofsTC = x;
      _yofsTC = y;
   }
   //-------------------------------------------------------------
   public void setMag(double mag) {
      _mag = mag;
   }

   //-------------------------------------------------------------
   public double getMag() {
      return _mag;
   }

   //-------------------------------------------------------------
   public Point getPos() {
      return _pos;
   }

   //-------------------------------------------------------------
   public int getGreyVal() {
      return _greyVal;
   }

   //-------------------------------------------------------------
   public Dimension getImgSize() {
      return _imgSize;
   }

   //-------------------------------------------------------------
   public BufferedImage getGreyBufferedImage() {
      return _bufImage;
   }

   public BufferedImage [] getColorBufImageArray(){
     return _abufImage;
   }

   public double[] getColorBufImgX() {
     return _axofs;
   }

   public double[] getColorBufImgY() {
     return _ayofs;
     }

   //-------------------------------------------------------------
   public BufferedImage getComponentBufferedImage(boolean showInterSecLines) {
     _showInterSecLines = showInterSecLines;
     getComponentBufferedImage();
   }
   //-------------------------------------------------------------
   public BufferedImage getComponentBufferedImage() {
     BufferedImage _compImage =
         new BufferedImage(_bufImage.getWidth(), _bufImage.getHeight(),
                           BufferedImage.TYPE_INT_RGB);
     Graphics2D g = _compImage.createGraphics();
     g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER));

     drawGreyImage(g);
     drawOverlay(g);
     drawAnatomy(g);
     if (_showInterSecLines) drawIntersection(g);
     return _compImage;
   }

   //-------------------------------------------------------------
   public WritableRaster getRaster() {
      return _ras;
   }

   //-------------------------------------------------------------
   // handle all _objects that are interested in changes
   //-------------------------------------------------------------
   // keep track of all the listeners to this 'model'
   protected EventListenerList changeListeners =
      new EventListenerList();

   //-------------------------------------------------------------
   // add a listener to the register
   public void addChangeListener(ChangeListener x) {
      changeListeners.add (ChangeListener.class, x);

      // bring it up to date with current state
      x.stateChanged(new ChangeEvent(this));
   }

   //-------------------------------------------------------------
   // remove a listener from the register
   public void removeChangeListener(ChangeListener x) {
      changeListeners.remove (ChangeListener.class, x);
   }

   //-------------------------------------------------------------
   private ChangeEvent ce;
   private Object[] listeners;
   private ChangeListener cl;
   public void fireChange() {
      // Create the event:
      ce = new ChangeEvent(this);
      // Get the listener list
      listeners = changeListeners.getListenerList();
      // Process the listeners last to first
      // List is in pairs, Class and instance
      for (int i
	    = listeners.length-2; i >= 0; i -= 2) {
	 if (listeners[i] == ChangeListener.class) {
	    ChangeListener cl = (ChangeListener)listeners[i+1];
	    cl.stateChanged(ce);
	 }
      }
   } // fireChange

} // class WlzImgView
