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

/**
 *   Responsible for drawing the screen image.
 *   <br>Uses RGB colour model with alpha (transparency).
 */
public class WlzImgView extends Component {

   /**   toggles the display of debugging messages */
   private final boolean _debug = false;

   /**  The underlying 2D Woolz object */
   private WlzObject _obj = null;

   /**
    *   Woolz data structure representing the bounding box
    *   of a 2D grey-level image
    *   contains xmin,ymin,xmax,ymax.
    */
   private WlzIBox2 _bBox = null;

   /**
    *   Origin of a 2D section through a grey-level Woolz object.
    *   It is Xmin,Ymin of the 2D bounding box and is
    *   not necessarily 0,0.
    */
   private WlzIVertex2 _org = null;

   /**
    *   Origin of a 2D section through a Woolz object,
    *   the result of a <em>mouse-click anatomy</em> operation,
    *   in the same space as the grey-level section.
    */
   private WlzIVertex2 _oorg = null;

   /**
    *   Origin of a 2D Woolz object,
    *   the result of a thresholding operation on a grey-level section,
    *   in the same space as the grey-level section.
    */
   private WlzIVertex2 _torg = null;

   /**
    *   Origin of a 2D section through a Woolz object,
    *   the result of an <em>anatomy menu</em> selection,
    *   in the same space as the grey-level section.
    */
   private WlzIVertex2 _aorg[] = null;

   /**
    *   Woolz data structure for a 2D grey-level Woolz object.
    *   (<em>grey value work space</em>).
    */
   private WlzGreyValueWSpace _gVWSp = null;

   /**
    *   Java2D data structure for a 2D grey-level image.
    *   Represents the grey-level image that is drawn to screen.
    */
   private BufferedImage _bufImage = null;

   /**
    *   Java2D data structure for an inverted 2D grey-level image.
    *   The same as _bufImage with all grey values inverted
    *   using: val = 255 - val.
    */
   private BufferedImage _invBufImage = null;

   /**
    *   Java2D data structure for a 2D monochrome image (magenta),
    *   representing the <em>mouse-click anatomy</em> image
    *   that is drawn to screen.
    */
   private BufferedImage _obufImage = null;

   /**
    *   Java2D data structure for a 2D monochrome image,
    *   representing the <em>threshold</em> image
    *   that is drawn to screen.
    */
   private BufferedImage _tbufImage = null;

   /**
    *   Java2D data structure for a 2D monochrome image,
    *   representing the <em>anatomy from menu</em> image
    *   that is drawn to screen.
    */
   private BufferedImage _abufImage[] = null;

   /**  Scale factor for the screen image */
   private double _mag;

   /**
    *   X coordinate of screen offset for
    *   <em>intersection lines</em>.
    */
   private double _ixofs;

   /**
    *   Y coordinate of screen offset for
    *   <em>intersection lines</em>.
    */
   private double _iyofs;

   /**
    *   X coordinate of screen offset for
    *   <em>mouse-click anatomy</em> image.
    */
   private double _oxofs;

   /**
    *   Y coordinate of screen offset for
    *   <em>mouse-click anatomy</em> image.
    */
   private double _oyofs;

   /**
    *   X coordinate of screen offset for
    *   <em>threshold</em> image.
    */
   private double _txofs;

   /**
    *   Y coordinate of screen offset for
    *   <em>threshold</em> image.
    */
   private double _tyofs;

   /**
    *   Array of X coordinates for screen offset of
    *   all of the <em>anatomy from menu</em> images.
    */
   private double _axofs[] = null;

   /**
    *   Array of Y coordinates for screen offset of
    *   all of the <em>anatomy from menu</em> images.
    */
   private double _ayofs[] = null;

   /**
    *   X coordinate of screen offset for
    *   <em>grey-level</em> image.
    */
   private double _xofsGL;

   /**
    *   Y coordinate of screen offset for
    *   <em>grey-level</em> image.
    */
   private double _yofsGL;

   /**
    *   X coordinate of screen offset for
    *   <em>threshold constraint</em>.
    */
   private double _xofsTC;

   /**
    *   Y coordinate of screen offset for
    *   <em>threshold constraint</em>.
    */
   private double _yofsTC;

   /**
    *   X coordinate of screen offset for
    *   <em>fixed point</em>.
    */
   private double _xofsFP;

   /**
    *   Y coordinate of screen offset for
    *   <em>fixed point</em>.
    */
   private double _yofsFP;

   /**
    *   X coordinate of screen offset for
    *   <em>fixed line</em>.
    */
   private double _xofsAP;

   /**
    *   Y coordinate of screen offset for
    *   <em>fixed line</em>.
    */
   private double _yofsAP;


   /**
    *   String representation of x,y coordinates and grey-level value
    *   at last mouse click / drag.
    *   This is in the same space as the grey-level section,
    *   (not necessarily the same as the screen coordinates).
    */
   private String _objStats;

   /**  Grey value (0 to 255) at last mouse click / drag */
   private int _greyVal;

   /**
    *   x,y coordinates at last mouse click / drag.
    *   This is in the same space as the grey-level section,
    *   (not necessarily the same as the screen coordinates).
    */
   private Point _pos = null;

   /**
    *   Java data structure representing the bounding box
    *   of a 2D <em>grey-level</em> image.
    */
   private Rectangle _bRec = null;

   /**
    *   Java data structure representing the bounding box
    *   of a 2D <em>mouse-click anatomy</em> image.
    */
   private Rectangle _obRec = null;

   /**
    *   Java data structure representing the bounding box
    *   of a <em>threshold</em> image.
    */
   private Rectangle _tbRec = null;

   /**
    *   Java data structure representing
    *   an array of bounding boxes for
    *   all <em>anatomy from menu</em> images.
    */
   private Rectangle _abRec[] = null;

   /**
    *   Java data structure representing
    *   the size of the image in pixels.
    */
   private Dimension _imgSize = null;

   /**
    *   Java2D data structure representing
    *   RGB colour space with transparency (alpha).
    */
   private DirectColorModel _dcm = null;

   /**
    *   Java2D data structure representing
    *   editable pixel data for
    *   the <em>grey-level</em> image.
    */
   private WritableRaster _ras = null;

   /**
    *   Java2D data structure representing
    *   editable pixel data for
    *   the <em>inverted grey-level</em> image.
    */
   private WritableRaster _invRas = null;

   /**
    *   Java2D data structure representing
    *   editable pixel data for
    *   the <em>mouse-click anatomy</em> image.
    */
   private WritableRaster _oras = null;

   /**
    *   Java2D data structure representing
    *   editable pixel data for
    *   the <em>threshold</em> image.
    */
   private WritableRaster _tras = null;

   /**
    *   Java2D data structure representing
    *   editable pixel data for
    *   all the <em>anatomy from menu</em> images.
    */
   private WritableRaster _aras[] = null;

   /**
    *   Java2D data structure containing an
    *   int data array used to construct the
    *   WritableRaster for <em>grey-level</em> images.
    */
   private DataBuffer _dataBuf = null;

   /**
    *   Java2D data structure containing an
    *   int data array used to construct the
    *   WritableRaster for <em>inverted grey-level</em> images.
    */
   private DataBuffer _invDataBuf = null;

   /**
    *   Java2D data structure containing an
    *   int data array used to construct the
    *   WritableRaster for <em>mouse-click anatomy</em> images.
    */
   private DataBuffer _odataBuf = null;

   /**
    *   Java2D data structure containing an
    *   int data array used to construct the
    *   WritableRaster for <em>threshold</em> images.
    */
   private DataBuffer _tdataBuf = null;

   /**
    *   Java2D data structure containing
    *   int data arrays used to construct the
    *   WritableRasters for all <em>anatomy from menu</em> images.
    */
   private DataBuffer _adataBuf[] = null;

   /**
    *   Java2D data structure containing a
    *   model which knows how to retrieve pixel samples
    *   from a <em>grey-level</em> DataBuffer.
    */
   private SinglePixelPackedSampleModel _sppsm = null;

   /**
    *   Java2D data structure containing a
    *   model which knows how to retrieve pixel samples
    *   from a <em>mouse-click anatomy</em> DataBuffer.
    */
   private SinglePixelPackedSampleModel _osppsm = null;

   /**
    *   Java2D data structure containing a
    *   model which knows how to retrieve pixel samples
    *   from a <em>threshold</em> DataBuffer.
    */
   private SinglePixelPackedSampleModel _tsppsm = null;

   /**
    *   Java2D data structure containing
    *   models which know how to retrieve pixel samples
    *   from all the <em>anatomy from menu</em> DataBuffers.
    */
   private SinglePixelPackedSampleModel _asppsm[] = null;

   /**
    *   The width of the colour model,
    *   8 bits each for R G B and alpha.
    */
   private int _nbits = 32;

   /**
    *   Array of bitmasks for the colour model.
    *   One each for R G B and alpha.
    *   These allow the 8 bit RGB and alpha values
    *   to be extracted from a 32 bit integer.
    */
   private int _masks[] = {0xff0000, 0xff00, 0xff, 0xff000000}; // (R G B alpha)

   /**
    *   Array containing data derived from a Woolz object for
    *   a <em>grey-level</em> image.
    */
   private int _imageData[] = null;

   /**
    *   Array containing data derived from a Woolz object for
    *   an <em>inverted grey-level</em> image.
    */
   private int _invImageData[] = null;

   /**
    *   Array containing data derived from a Woolz object for
    *   a <em>mouse-click anatomy</em> image.
    */
   private int _oimageData[] = null;

   /**
    *   Array containing data derived from a Woolz object for
    *   a <em>threshold</em> image.
    */
   private int _timageData[] = null;

   /**
    *   Array containing data derived from Woolz objects for
    *   all <em>anatomy from menu</em> images.
    */
   private int _aimageData[][] = null;


   /**
    *   Represents a polygon used to constrain thresholding.
    */
   private GeneralPath _threshConstraintPath = null;

   /**
    *   Number of points in a <em>threshold constraint</em> polygon.
    */
   private int _npts = 0;

   /**
    *   Local copy of the number of rows in AnatKey.
    */
   private int _num;

   /**
    *   Radius of the circle which denotes the <em>fixed point</em>.
    */
   private int _fpr = 2; // fixed point indicator radius

   /**
    *   Radius of the circle which denotes the <em>2nd fixed point</em>.
    *   Obsolete.
    */
   private int _apr = 2; // axis point indicator radius


   /**
    *   Array of booleans indicating the visibility
    *   of each of the <em>anatomy from menu</em> components.
    */
   private boolean _viz[]; // visibility of anatomy objects

   /**
    *   Collection of <em>Line2D.Double</em> elements representing
    *   the intersection of other SectionViewers with this one.
    */
   private Vector _intersectionVec = null;

   /**
    *   Collection of <em>Color</em> elements representing the colour of
    *   the corresponding intersection lines in <em>_intersectionVec</em>.
    */
   private Vector _interColVec = null;

   /**
    *   Collection of <em>Line2D.Double</em> elements representing
    *   a cross centred on the <em>fixed point</em>.
    *   A cross pattern is no longer used to indicate the fixed point,
    *   the first element is used to obtain the centre point for a circle.
    */
   private Vector _fixedPointVec = null;

   /**
    *   Collection of <em>Line2D.Double</em> elements representing
    *   a cross centred on the <em>2nd fixed point</em>.
    *   Obsolete.
    */
   private Vector _axisPointVec = null;

   /**
    *   x,y, coordinates of the <em>2nd fixed point</em>.
    */
   private double _axisPointArr[] = null;

   /**  True if <em>fixed point</em> is to be drawn. */
   private boolean _fixedPoint;

   /**  True if <em>2nd fixed point</em> is to be drawn. */
   private boolean _axisPoint;

   /**  True if <em>fixed line</em> is to be drawn. */
   private boolean _axis;

   /**  True if <em>intersection lines</em> are to be drawn. */
   private boolean _intersection;

   /**  True if <em>anatomy from menu</em> components are to be drawn. */
   private boolean _anatomy;

   /**  True if <em>thresholded</em> region is to be drawn. */
   private boolean _threshold;

   /**  True if <em>threshold constraint</em> polygon is to be drawn. */
   private boolean _threshConstraint;

   /**  True if <em>mouse-click anatomy</em> is to be drawn. */
   private boolean _overlay;

   /**  The colour of the mouse-click anatomy overlay */
   private Color _mcCol;

   /**
    *   True if <em>grey-level</em> image
    *   is to be drawn with inverted grey values.
    */
   private boolean _inverted = false;

//...................................
   /**  Required by Tie Point application. */
   private boolean _showInterSecLines = false;

   /**  Required by Tie Point application. */
   private boolean _tiepoint = true;

   /**  Required by Tie Point application. */
   private Vector tps = null;

   /**  Required by Tie Point application. */
   private Vector tpsCol = null;

   //-------------------------------------------------------------
   // constructor
   /**
    *   Constructs a new WlzImgView.
    */
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

      _num = 0;

      _intersectionVec = new Vector();
      _interColVec = new Vector();
      _interColVec.add(Color.red);
      _fixedPointVec = new Vector();
      _axisPointVec = new Vector();

      _mcCol = new Color(255,0,255,125); // default is magenta semi-transparent
   }

   //-------------------------------------------------------------
   /**
    *   Generates a <em>32 bit single pixel packed int</em> array
    *   from the given 2D grey-level Woolz object,
    *   and creates a Java2D DirectColorModel that knows
    *   how to extract RGB and alpha information from
    *   a 32 bit single pixel packed integer.
    *   <ul>
    *   <li>Woolz data is converted to a 1 dimensional
    *   unsigned byte (8 bit) array whose length =
    *   Width x Length of the 2D bounding box.</li>
    *   <li>Each unsigned byte is converted to a 32 bit 2s complement int.</li>
    *   <li>The 8-bit pattern from the Woolz data is bit shifted by 8 and 16
    *   and combined to create a grey value.</li>
    *   <li>the alpha value (255) is shifted by 24 and combined with
    *   the RGB to create a 32 bit single pixel packed integer.</li>
    *   </ul>
    *   @param newObj the 2D grey-level Woolz object.
    */
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
      //setSize(_imgSize);


      _dcm = new DirectColorModel(
	    _nbits,
	    _masks[0],
	    _masks[1],
	    _masks[2],
	    _masks[3]);
      //System.out.println(_dcm.toString());

      _imageData = new int[_bRec.width * _bRec.height];
      _invImageData = new int[_bRec.width * _bRec.height];
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
            // make the inverted data
            ii = 255 - ii;
            _invImageData[(idY * _bRec.width) + idX] =
               ii | (ii << 8) | (ii << 16) | (255 << 24);
	 } // for
      } // for

      makeBufImage();
   }

   //-------------------------------------------------------------
   /**
    *   Generates a Java2D BufferedImage from the
    *   <em>32 bit single pixel packed int</em> array
    *   created by setWlzObj().
    *   <ul>
    *   <li>Creates a Java2D DataBuffer from the array and
    *   the 2D bounding box dimensions.</li>
    *   <li>Creates a Java2D SinglePixelPackedSampleModel
    *   which describes how to retrieve data from the DataBuffer.</li>
    *   <li>Creates a Java2D WritableRaster given the above
    *   SinglePixelPackedSampleModel and the DataBuffer.</li>
    *   <li>Creates a BufferedImage using the DirectColorModel,
    *   from setWlzObj(), and the above WritableRaster.</li>
    *   </ul>
    */
   public void makeBufImage() {

      _dataBuf = new DataBufferInt(_imageData,
	    _bRec.width * _bRec.height);
      if(_dataBuf == null) {
	 System.out.println("_dataBuf is null");
      }
      // make the inverted one
      _invDataBuf = new DataBufferInt(_invImageData,
            _bRec.width * _bRec.height);

      _sppsm = new SinglePixelPackedSampleModel(
	    DataBuffer.TYPE_INT,
	    _bRec.width,
	    _bRec.height,
	    _masks);

      _ras = Raster.createWritableRaster(_sppsm,
	    _dataBuf,
	    null);
      // make the inverted one
      _invRas = Raster.createWritableRaster(_sppsm,
            _invDataBuf,
            null);

      try {
	 _bufImage = new BufferedImage(_dcm, _ras, false, null);
	 _invBufImage = new BufferedImage(_dcm, _invRas, false, null);
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
   /**
    *   Generates a <em>32 bit single pixel packed int</em> array
    *   from the given <em>mouse-click anatomy</em> 2D Woolz object.
    *   <ul>
    *   <li>Woolz data is converted to a 1 dimensional
    *   unsigned byte (8 bit) array whose length =
    *   Width x Length of the anatomy object's bounding box.</li>
    *   <li>Each unsigned byte is converted to a 32 bit 2s complement int.</li>
    *   <li>Black pixels (background) are converted to
    *   fully transparent white.</li>
    *   <li>Other pixels are converted to semi-transparent Magenta.</li>
    *   </ul>
    *   @param newObj the 2D Woolz object.
    */
   public void setOverlayObj(WlzObject newObj) throws WlzException {
      WlzObject obj = newObj;
      WlzIBox2 bBox = null;
      Dimension imgSize = new Dimension();
      //int imageData[] = null;
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
		     (_mcCol.getAlpha() << 24) |
		     (_mcCol.getRed() << 16) |
		     (_mcCol.getGreen() << 8) |
		     (_mcCol.getBlue());
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
   /**
    *   Generates a Java2D BufferedImage from the
    *   <em>32 bit single pixel packed int</em> array
    *   created by setOverlayObj().
    *   <ul>
    *   <li>Creates a Java2D DataBuffer from the array and
    *   the anatomy object's bounding box dimensions.</li>
    *   <li>Creates a Java2D SinglePixelPackedSampleModel
    *   which describes how to retrieve data from the DataBuffer.</li>
    *   <li>Creates a Java2D WritableRaster given the above
    *   SinglePixelPackedSampleModel and the DataBuffer.</li>
    *   <li>Creates a BufferedImage using the DirectColorModel,
    *   from setWlzObj(), and the above WritableRaster.</li>
    *   </ul>
    */
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
   /**
    *   Generates a <em>32 bit single pixel packed int</em> array
    *   from the given <em>threshold</em> Woolz object.
    *   <ul>
    *   <li>Woolz data is converted to a 1 dimensional
    *   unsigned byte (8 bit) array whose length =
    *   Width x Length of the threshold object's bounding box.</li>
    *   <li>Each unsigned byte is converted to a 32 bit 2s complement int.</li>
    *   <li>White pixels (background) are made fully transparent.</li>
    *   <li>Other pixels are converted to semi-transparent Green.</li>
    *   </ul>
    *   @param newObj the Woolz object.
    */
   public void setThresholdObj(WlzObject newObj) throws WlzException {
      WlzObject obj = newObj;
      WlzIBox2 bBox = null;
      Dimension imgSize = new Dimension();
      //int imageData[] = null;
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
   /**
    *   Generates a Java2D BufferedImage from the
    *   <em>32 bit single pixel packed int</em> array
    *   created by setThresholdObj().
    *   <ul>
    *   <li>Creates a Java2D DataBuffer from the array and
    *   the threshold object's bounding box dimensions.</li>
    *   <li>Creates a Java2D SinglePixelPackedSampleModel
    *   which describes how to retrieve data from the DataBuffer.</li>
    *   <li>Creates a Java2D WritableRaster given the above
    *   SinglePixelPackedSampleModel and the DataBuffer.</li>
    *   <li>Creates a BufferedImage using the DirectColorModel,
    *   from setWlzObj(), and the above WritableRaster.</li>
    *   </ul>
    */
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
   /**
    *   Generates a collection of
    *   <em>32 bit single pixel packed int</em> arrays
    *   from the given <em>anatomy from menu</em> 2D Woolz objects.
    *   <ul>
    *   <li>Woolz data is converted to a 1 dimensional
    *   unsigned byte (8 bit) array whose length =
    *   Width x Length of the anatomy object's bounding box.</li>
    *   <li>Each unsigned byte is converted to a 32 bit 2s complement int.</li>
    *   <li>Black pixels (background) are converted to
    *   fully transparent white.</li>
    *   <li>Other pixels are converted to a semi-transparent colour.
    *   obtained from the Anatomy Key.</li>
    *   </ul>
    *   @param objs the array of 2D Woolz objects.
    *   @param viz an array of booleans indicating the visibility
    *   of each of the anatomy components.
    */
   public void setAnatomyObj(WlzObject[] objs,
                             boolean[] viz,
			     Color[] col)
      throws WlzException {

	 _viz = viz;

	 if(objs.length != _num) {
	    _num = objs.length;
	    initArrays();
	 }
	 for(int i = 0; i < _num; i++) {
	    if(_viz[i] == false) {
	       _aorg[i] = null;
	       _abufImage[i] = null;
	       continue;
	    }
	    WlzObject obj = objs[i];
	    if(obj != null) {
	       WlzIBox2 bBox = null;
	       Dimension imgSize = new Dimension();
	       //int imageData[] = null;
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
                           ii = col[i].getBlue() |
                                col[i].getGreen() << 8 |
                                col[i].getRed() << 16 |
                                125 << 24;

			   _aimageData[i][(idY * _abRec[i].width) + idX] = ii;
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
   /**
    *   Generates a Java2D BufferedImage from a
    *   <em>32 bit single pixel packed int</em> array
    *   created by setAnatomyObj().
    *   <ul>
    *   <li>Creates a Java2D DataBuffer from the array and
    *   the anatomy object's bounding box dimensions.</li>
    *   <li>Creates a Java2D SinglePixelPackedSampleModel
    *   which describes how to retrieve data from the DataBuffer.</li>
    *   <li>Creates a Java2D WritableRaster given the above
    *   SinglePixelPackedSampleModel and the DataBuffer.</li>
    *   <li>Creates a BufferedImage using the DirectColorModel,
    *   from setWlzObj(), and the above WritableRaster.</li>
    *   </ul>
    *   @param indx index into the collection of
    *   <em>32 bit single pixel packed int</em> arrays.
    */
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
   /**
    *   Creates and initialises arrays required for
    *   anatomy from menu.
    */
   private void initArrays() {
      _abufImage = null;
      _adataBuf = null;
      _asppsm = null;
      _aimageData = null;
      _aorg = null;
      _axofs = null;
      _ayofs = null;
      _abRec = null;
      _aras = null;

      _abufImage = new BufferedImage[_num];
      _adataBuf = new DataBuffer[_num];
      _asppsm = new SinglePixelPackedSampleModel[_num];
      _aimageData = new int[_num][];
      _aorg = new WlzIVertex2[_num];
      _axofs = new double[_num];
      _ayofs = new double[_num];
      _abRec = new Rectangle[_num];
      _aras = new WritableRaster[_num];

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
   }

   //-------------------------------------------------------------
   // makes a GeneralPath from Vector of Floats
   /**
    *   Creates a Java2d GeneralPath representing
    *   a polygonal region in which to threshold.
    *   @param xPoints Collection of x coordinates of points
    *   defining the polygonal constraint region.
    *   @param yPoints Collection of y coordinates of points
    *   defining the polygonal constraint region.
    */
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
   /**
    *   Makes the threshold constraint a closed polygon by
    *   joining the first and last points.
    */
   public void closeThreshConstraint() {
      _threshConstraintPath.closePath();
   }
   //-------------------------------------------------------------
   /**
    *   Draws a composite image to a graphics device,
    *   normally the screen.
    *   The composite image may comprise
    *   <ul>
    *   <li>grey-level section</li>
    *   <li>mouse-click anatomy</li>
    *   <li>anatomy from menu</li>
    *   <li>intersection lines</li>
    *   <li>fixed point</li>
    *   <li>fixed line</li>
    *   <li>threshold constraint</li>
    *   <li>threshold region</li>
    *   </ul>
    *   If the SectionViewer is being used in a Tie Point application
    *   the composite image may also contain tie points.
    *   @param g Java Graphics object which encapsulates
    *   state information needed for the basic rendering operations
    *   that Java supports.
    */
   public void paint(Graphics g) {
      Graphics2D g2 = (Graphics2D)g;
      g2.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER));
      g2.scale(_mag, _mag);
      drawGreyImage(g2);
      drawOverlay(g2);
      drawAnatomy(g2);
      drawIntersection(g2);
      drawFixedPoint(g2);
      //drawAxisPoint(g2);
      drawAxis(g2);
      // temporarily disabled ... don't remove
      drawThreshold(g2);
      drawThreshConstraint(g2);
      //--------------------------------------
      drawTiePoint(g2);
      drawWrapedTiePoint(g2);
   }

   //-------------------------------------------------------------
   /**
    *   Draws the <em>grey-level section</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
   public void drawGreyImage(Graphics2D g2) {
      if(_bufImage != null) {
	 g2.translate(_xofsGL, _yofsGL);
         if(_inverted) {
            g2.drawImage(_invBufImage, 0, 0, null);
         } else {
            g2.drawImage(_bufImage, 0, 0, null);
         }
	 g2.translate(-_xofsGL, -_yofsGL);
      }
   }

   //-------------------------------------------------------------
   /**
    *   Draws the <em>mouse-click anatomy</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
   public void drawOverlay(Graphics2D g2) {
      if (_overlay == false) return;
      if(_obufImage != null) {
	 g2.translate(_oxofs, _oyofs);
	 g2.drawImage(_obufImage, 0, 0, null);
	 g2.translate(-_oxofs, -_oyofs);
      }
   }

   //-------------------------------------------------------------
   /**
    *   Draws the <em>thresholded region</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
   public void drawThreshold(Graphics2D g2) {
      if (_threshold == false) return;
      if(_tbufImage != null) {
	 g2.translate(_txofs, _tyofs);
	 g2.drawImage(_tbufImage, 0, 0, null);
	 g2.translate(-_txofs, -_tyofs);
      }
   }

   //-------------------------------------------------------------
   /**
    *   Draws the <em>anatomy from menu components</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
   public void drawAnatomy(Graphics2D g2) {
      if (_anatomy == false) return;

      for(int i=0; i<_num; i++) {
	 if(_viz[i] == false) continue;
	 if(_abufImage[i] != null) {
	    g2.translate(_axofs[i], _ayofs[i]);
	    g2.drawImage(_abufImage[i], 0, 0, null);
	    g2.translate(-_axofs[i], -_ayofs[i]);
	 }
      }
   }

   //-------------------------------------------------------------
   /**
    *   Draws the <em>intersection lines</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
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
   /**
    *   Draws the <em>fixed point indicator</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
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
      line = (Line2D.Double)_fixedPointVec.elementAt(0);
      int x = (int)line.getX1() - _fpr;
      int y = (int)line.getY1() - _fpr;
      g2.drawOval(x,y,2*_fpr,2*_fpr);
/*
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
*/
      g2.translate(-_xofsFP, -_yofsFP);
      g2.setColor(colOrig);
   }

   //-------------------------------------------------------------
   /**
    *   Draws the <em>2nd fixed point indicator</em> to a graphics device,
    *   normally the screen.
    *   Obsolete.
    *   @param g Java2D Graphics object.
    */
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
   /**
    *   Draws the <em>fixed line</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
   public void drawAxis(Graphics2D g2) {

      // draw a line between fixed point and axis point
      if(_axis == false) return;
      if((_fixedPointVec == null) || (_fixedPointVec.size() == 0)) return;
      if(_axisPointArr == null) return;

      Color colOrig = g2.getColor();
      Line2D.Double line = null;

      g2.setColor(Color.cyan);
      g2.translate(_xofsAP, _yofsAP);

      int x1 = (int)_axisPointArr[0];
      int y1 = (int)_axisPointArr[1];

      line = (Line2D.Double)_fixedPointVec.elementAt(0);
      int x2 = (int)line.getX1();
      int y2 = (int)line.getY1();

      g2.drawLine(x1,y1,x2,y2);

      g2.translate(-_xofsAP, -_yofsAP);
      g2.setColor(colOrig);
   }

   //-------------------------------------------------------------
   /**
    *   Draws the <em>threshold constraint polygon</em> to a graphics device,
    *   normally the screen.
    *   @param g Java2D Graphics object.
    */
   public void drawThreshConstraint(Graphics2D g2) {
      if(_threshConstraint && (_threshConstraintPath != null)) {
	 g2.setColor(Color.red);
	 g2.draw(_threshConstraintPath);
      }
   }

   //-------------------------------------------------------------
   /**
    *   Removes <em>mouse-click anatomy</em> from the composite image.
    */
   public void clearOverlay() {
      _obufImage = null;
      _overlay = false;
   }

   //-------------------------------------------------------------
   /**
    *   Toggles display of the <em>fixed point</em> representation.
    */
   public void enableFixedPoint(boolean state) {
      _fixedPoint = state;
   }

   //-------------------------------------------------------------
   /**
    *   Toggles display of the <em>2nd fixed point</em> representation.
    *   Obsolete.
    */
   public void enableAxisPoint(boolean state) {
      _axisPoint = state;
   }

   //-------------------------------------------------------------
   /**
    *   Toggles display of the <em>fixed line</em>.
    */
   public void enableAxis(boolean state) {
      _axis = state;
   }

   //-------------------------------------------------------------
   /**
    *   Removes <em>mouse-click anatomy</em> from the composite image.
    */
   public void clearIntersection() {
      _intersection = false;
   }

   //-------------------------------------------------------------
   /**
    *   Removes <em>mouse-click anatomy</em> from the composite image.
    */
   public void clearThreshold() {
      _tbufImage = null;
      _threshold = false;
   }

   //-------------------------------------------------------------
   /**
    *   Toggles display of the <em>threshold constraint</em>.
    */
   public void enableThreshConstraint(boolean state) {
      _threshConstraint = state;
   }

   //-------------------------------------------------------------
   /**
    *   Removes <em>threshold constraint</em> from the composite image.
    */
   public void clearThreshConstraint() {
      _threshConstraintPath = null;
   }

   //-------------------------------------------------------------
   /**
    *   Removes <em>anatomy from menu</em> from the composite image.
    */
   public void clearAnatomy() {
      for(int i=0; i<_num; i++) {
	 _abufImage[i] = null;
      }
      _anatomy = false;
   }

   //-------------------------------------------------------------
   /**
    *   Returns the GeneralPath representing a threshold constraint polygon.
    */
   public GeneralPath getThreshConstraint() {
      return _threshConstraintPath;
   }

   //-------------------------------------------------------------
   /**
    *   Stores grey-level value and  position at the last
    *   significant mouse event in the coordinate space of
    *   a 2D grey-level section.
    *   (Not necessarily the screen coordinate system.)
    *   Significant mouse events are
    *   <ul>
    *   <li>press</li>
    *   <li>drag</li>
    *   <li>release</li>
    *   <li>click</li>
    *   </ul>
    *   @param relX X coordinate of the screen position of the mouse event.
    *   @param relY Y coordinate of the screen position of the mouse event.
    */
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
   /**
    *   Returns the offset of a <em>grey-level section</em>
    *   from 0,0 (the screen origin).
    *   The origin of a <em>grey-level section</em> is
    *   obtained from its Woolz bounding box.
    *   @return  the offset as a Java Point.
    */
   public Point getOffSet(){
     return new Point(_bBox.xMin, _bBox.yMin);
   }

   //-------------------------------------------------------------
   /**
    *   Causes a <em>grey-level section</em> to be drawn
    *   with its grey-level values inverted.
    *   i.e. newGreyVal = 255 - greyVal.
    *   @param state true if inversion is to occur.
    */
   public void setInverted(boolean state) {
      _inverted = state;
   }

   //-------------------------------------------------------------
   /**
    *   Copies an array of Java Line2D.Double objects to a local Vector.
    *   @param lines the array of Line2D.Double objects.
    */
   public void setIntersectionVec(Line2D.Double[] lines) {

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
   /**
    *   Copies an array of Java Color objects to a local Vector.
    *   @param lines the array of Color objects.
    */
   public void setInterColVec(Color[] cols) {

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
   /**
    *   Required by Tie Point application.
    */
   void drawTiePoint(Graphics2D g){
     if (!_tiepoint) return;
     if (null == tps) return;

     Color orgColor = g.getColor();
     g.translate(_xofsFP, _yofsFP);
     g.scale(1/_mag, 1/_mag);
     for (int i = 0; i < tps.size(); i++){
       Point p = (Point) tps.get(i);
       g.setColor((Color)tpsCol.get(i));
       g.fillOval((int)(_mag*p.x) - 4, (int)(_mag*p.y) - 4, 8, 8);
     }
     g.scale(_mag, _mag);
     g.setColor(orgColor);
     g.translate(-_xofsFP, -_yofsFP);
   }

   private Vector itsP = null;
   void drawWrapedTiePoint(Graphics2D g){
     if (null == itsP) return;

     Color orgColor = g.getColor();
     g.translate(_xofsFP, _yofsFP);
     g.scale(1/_mag, 1/_mag);
     g.setColor(Color.orange);
     Point[] p = new Point[itsP.size()];
     for (int i = 0; i < p.length; i++){
       p[i] = (Point) itsP.get(i);
       g.fillOval((int)(_mag*p[i].x) - 1, (int)(_mag*p[i].y) - 1, 2, 2);
     }
     for (int i = 0; i < p.length; i+=2){
       g.drawLine((int)(_mag*p[i].x),(int)(_mag*p[i].y),
                  (int)(_mag*p[i+1].x),(int)(_mag*p[i+1].y));
     }

     g.scale(_mag, _mag);
     g.setColor(orgColor);
     g.translate(-_xofsFP, -_yofsFP);
   }

   //-------------------------------------------------------------
   /**
    *   Required by Tie Point application.
    */
   public void setTiePoint(Vector tps, Vector tpsCol){
     this.tps = tps;
     this.tpsCol = tpsCol;
     setFixedPointOffsets();
   }


   public void setWrapedTiePoint(Vector itsP){
     this.itsP = itsP;
     setFixedPointOffsets();
     this.repaint();
   }
   //-------------------------------------------------------------
   /**
    *   Required by Tie Point application.
    */
   public void enableTiePoint(boolean state) {
         _tiepoint = state;
   }
   //-------------------------------------------------------------
   /**
    *   Generates a local collection of Java Line2D.Double objects
    *   to indicate <em>fixed point</em> position as a cross.
    *   Currently the fixed point is indicated by a circle.
    *   @param fpa the x,y coordinates of the fixed point.
    */
   public void setFixedPointVec(double[] fpa) {

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
   /**
    *   Copies the coordinates of the <em>2nd fixed point</em>
    *   to a local Vector.
    *   @param apa the x,y coordinates of the 2nd fixed point.
    */
   public void setAxisPointArr(double[] apa) {

      if(apa == null) return;

      if(_axisPointArr != null) {
         _axisPointArr = null;
      }
      _axisPointArr = new double[3];

      _axisPointArr[0] = apa[0];
      _axisPointArr[1] = apa[1];
      _axisPointArr[2] = apa[2]; // not used

      setAxisPointOffsets();
      _axis = true;
   }


   //-------------------------------------------------------------
   /**
    *   Generates a local collection of Java Line2D.Double objects
    *   to indicate <em>2nd fixed point</em> position as a cross.
    *   Currently the 2nd fixed point is not displayed.
    *   @param fpa the x,y coordinates of the 2nd fixed point.
    */
   public void setAxisPointVec(double[] apa) {

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
   /**
    *   Debugging method, prints the x,y coordinates of each end
    *   of a line.
    *   @param line the line to print.
    */
   public void printLine(Line2D.Double line) {
      System.out.println(Double.toString(line.getX1())+", "+
                         Double.toString(line.getY1())+", "+
			 Double.toString(line.getX2())+", "+
			 Double.toString(line.getY2()));
   }

   //-------------------------------------------------------------
   /**
    *   Sets the offset between the origin of the graphics display (screen)
    *   and the origin of the BufferedImage representing a
    *   <em>grey-level section</em>.
    *   Note: this is not the same as the offset between the
    *   <em>grey-level section's</em> bounding box and the screen because the
    *   raster representing the <em>grey-level section</em> is filled from
    *   0,0.
    */
   public void setGreyLevelOffsets() {

      _xofsGL = 0;
      _yofsGL = 0;

   }
   //-------------------------------------------------------------
   /**
    *   Sets the offset between the origin of the graphics display (screen)
    *   and the origin of the coordinate space for <em>intersection lines</em>.
    *   Should be 0,0.
    */
   public void setIntersectionOffsets() {
      _ixofs = 0;
      _iyofs = 0;
   }
   //-------------------------------------------------------------
   /**
    *   Sets the offset between the origin of the graphics display (screen)
    *   and the <em>fixed point</em>.
    *   The <em>fixed point</em> is given in the coordinate space of a
    *   <em>grey-level section</em> and therefore needs to be offset
    *   when displayed on screen.
    */
   public void setFixedPointOffsets() {
      double X = Math.abs((_bBox.xMax - _bBox.xMin));
      double Y = Math.abs((_bBox.yMax - _bBox.yMin));

      _xofsFP = X - _bBox.xMax;
      _yofsFP = Y - _bBox.yMax;
   }
   //-------------------------------------------------------------
   /**
    *   Returns the offset between the origin of the graphics display (screen)
    *   and the <em>fixed point</em>.
    *   The <em>fixed point</em> is given in the coordinate space of a
    *   <em>grey-level section</em> and therefore needs to be offset
    *   when displayed on screen.
    */
   public double[] getFixedPointOffsets() {
      double ret[] = new double[2];

      ret[0] = _xofsFP;
      ret[1] = _yofsFP;

      return ret;
   }
   //-------------------------------------------------------------
   /**
    *   Sets the offset between the origin of the graphics display (screen)
    *   and the <em>2nd fixed point</em>.
    *   The <em>2nd fixed point</em> is given in the coordinate space of a
    *   <em>grey-level section</em> and therefore needs to be offset
    *   when displayed on screen.
    */
   public void setAxisPointOffsets() {
      double X = Math.abs((_bBox.xMax - _bBox.xMin));
      double Y = Math.abs((_bBox.yMax - _bBox.yMin));

      _xofsAP = X - _bBox.xMax;
      _yofsAP = Y - _bBox.yMax;
   }
   //-------------------------------------------------------------
   /**
    *   Sets the offset between the origin of the graphics display (screen)
    *   and <em>mouse-click anatomy</em>.
    *   <em>Mouse-click anatomy</em> is in the coordinate space of a
    *   <em>grey-level section</em> and therefore needs to be offset
    *   when displayed on screen.
    */
   public void setOverlayOffsets() {
      _oxofs = Math.abs(_oorg.vtX - _org.vtX);
      _oyofs = Math.abs(_oorg.vtY - _org.vtY);
   }
   //-------------------------------------------------------------
   /**
    *   Sets the offset between the origin of the graphics display (screen)
    *   and a <em>threshold</em> region.
    *   A <em>threshold</em> region is in the coordinate space of a
    *   <em>grey-level section</em> and therefore needs to be offset
    *   when displayed on screen.
    */
   public void setThresholdOffsets() {
      _txofs = Math.abs(_torg.vtX - _org.vtX);
      _tyofs = Math.abs(_torg.vtY - _org.vtY);
   }
   //-------------------------------------------------------------
   /**
    *   Sets offsets between the origin of the graphics display (screen)
    *   and all the <em>anatomy from menu</em> components.
    *   <em>Anatomy from menu</em> regions are in the coordinate space of a
    *   <em>grey-level section</em> and therefore needs to be offset
    *   when displayed on screen.
    */
   public void setAnatomyOffsets() {

      for(int i=0; i<_num; i++) {
	 if(_aorg[i] != null) {
	    _axofs[i] = Math.abs(_aorg[i].vtX - _org.vtX);
	    _ayofs[i] = Math.abs(_aorg[i].vtY - _org.vtY);
	 }
      }
   }
   //-------------------------------------------------------------
   /**
    *   Sets offsets between the origin of the graphics display (screen)
    *   and a <em>threshold constraint</em> polygon.
    *   Should be 0,0.
    */
   public void setThreshConstraintOffsets() {
      _xofsTC = 0;
      _yofsTC = 0;
   }
   //-------------------------------------------------------------
   /**
    *   Sets the scale factor for the screen image.
    *   @param mag the scale factor.
    */
   public void setMag(double mag) {
      _mag = mag;
   }

   //-------------------------------------------------------------
   /**
    *   Gets the scale factor for the screen image.
    *   @return the scale factor.
    */
   public double getMag() {
      return _mag;
   }

   //-------------------------------------------------------------
   /**
    *   Returns the position of the last
    *   significant mouse event in the coordinate space of
    *   a 2D grey-level section.
    *   (Not necessarily the screen coordinate system.)
    *   Significant mouse events are
    *   <ul>
    *   <li>press</li>
    *   <li>drag</li>
    *   <li>release</li>
    *   <li>click</li>
    *   </ul>
    *   @return the position.
    */
   public Point getPos() {
      return _pos;
   }

   //-------------------------------------------------------------
   /**
    *   Returns the grey-level value at the last
    *   significant mouse event in the coordinate space of
    *   a 2D grey-level section.
    *   (Not necessarily the screen coordinate system.)
    *   Significant mouse events are
    *   <ul>
    *   <li>press</li>
    *   <li>drag</li>
    *   <li>release</li>
    *   <li>click</li>
    *   </ul>
    *   @return the grey-level value.
    */
   public int getGreyVal() {
      return _greyVal;
   }
   public void setShowInterSecLines(boolean flag) {
       _showInterSecLines = flag;
     }

   //-------------------------------------------------------------
   /**
    *   Returns the size of the image in pixels.
    *   @return the size.
    */
   public Dimension getImgSize() {
      return _imgSize;
   }

   //-------------------------------------------------------------
   /**
    *   Returns a BufferedImage representing a
    *   <em>grey-level section</em>.
    *   @return the BufferedImage.
    */
   public BufferedImage getGreyBufferedImage() {
      return _bufImage;
   }

   //-------------------------------------------------------------
   /**
    *   Returns an array of BufferedImages representing
    *   <em>anatomy from menu</em> components.
    *   @return the BufferedImage.
    */
   public BufferedImage [] getColorBufImageArray(){
     return _abufImage;
   }

   //-------------------------------------------------------------
   /**
    *   Returns an array representing the x offsets of
    *   <em>anatomy from menu</em> components.
    *   @return an array of x offsets.
    */
   public double[] getColorBufImgX() {
     return _axofs;
   }

   //-------------------------------------------------------------
   /**
    *   Returns an array representing the y offset of
    *   <em>anatomy from menu</em> components.
    *   @return an array of y offsets.
    */
   public double[] getColorBufImgY() {
     return _ayofs;
   }

   //-------------------------------------------------------------
   /**
    *   Returns a BufferedImage representing a composite image.
    *   The composite image may comprise
    *   <ul>
    *   <li>grey-level section</li>
    *   <li>mouse-click anatomy</li>
    *   <li>anatomy from menu</li>
    *   <li>intersection lines</li>
    *   <li>fixed point</li>
    *   <li>fixed line</li>
    *   <li>threshold constraint</li>
    *   <li>threshold region</li>
    *   </ul>
    *   If the SectionViewer is being used in a Tie Point application
    *   the composite image may also contain tie points.
    *   @return the BufferedImage.
    */
   public BufferedImage getComponentBufferedImage() {
     return getComponentBufferedImage(_showInterSecLines);
   }
   //-------------------------------------------------------------
   /**
    *   Returns a BufferedImage representing a composite image.
    *   The composite image may comprise
    *   <ul>
    *   <li>grey-level section</li>
    *   <li>mouse-click anatomy</li>
    *   <li>anatomy from menu</li>
    *   <li>intersection lines</li>
    *   <li>fixed point</li>
    *   <li>fixed line</li>
    *   <li>threshold constraint</li>
    *   <li>threshold region</li>
    *   </ul>
    *   If the SectionViewer is being used in a Tie Point application
    *   the composite image may also contain tie points.
    *   @param showInterSecLines an application such as JWlzViewer may want to
    *   display intersection lines in 3D even if they are not displayed in the
    *   2D composite image and vice-versa.
    *   @return the BufferedImage.
    */
   public BufferedImage getComponentBufferedImage(boolean showInterSecLines) {
     BufferedImage _compImage =
         new BufferedImage(_bufImage.getWidth(), _bufImage.getHeight(),
                           BufferedImage.TYPE_INT_RGB);
     Graphics2D g = _compImage.createGraphics();
     g.setComposite(AlphaComposite.getInstance(AlphaComposite.SRC_OVER));

     drawGreyImage(g);
     drawOverlay(g);
     drawAnatomy(g);
     if (showInterSecLines) drawIntersection(g);
     return _compImage;
   }

   //-------------------------------------------------------------
   /**
    *   Returns the Java2D WritableRaster that represents a
    *   <em>grey-level section</em>.
    *   @return the WritableRaster.
    */
   public WritableRaster getRaster() {
      return _ras;
   }

   //-------------------------------------------------------------
   /**
    *   Set the colour of mouse-click anatomy
    *   @param col, the colour of the overlay.
    */
   public void setMCCol(Color col) {
      _mcCol = col;
   }

   //-------------------------------------------------------------
   // handle all _objects that are interested in changes
   //-------------------------------------------------------------
   // keep track of all the listeners to this 'model'
   /**
    *   A list of ChangeListeners which are
    *   listening for events fired from the WlzImgView.
    */
   protected EventListenerList changeListeners =
      new EventListenerList();

   //-------------------------------------------------------------
   // add a listener to the register
   /**
    *   Adds a ChangeListener to the EventListenerList.
    */
   public void addChangeListener(ChangeListener x) {
      changeListeners.add (ChangeListener.class, x);

      // bring it up to date with current state
      x.stateChanged(new ChangeEvent(this));
   }

   //-------------------------------------------------------------
   // remove a listener from the register
   /**
    *   Removes a ChangeListener from the EventListenerList.
    */
   public void removeChangeListener(ChangeListener x) {
      changeListeners.remove (ChangeListener.class, x);
   }

   //-------------------------------------------------------------
   /**   An event that will be fired from WlzImgView */
   private ChangeEvent ce;

   /**  A local copy of the list of ChangeListeners */
   private Object[] listeners;

   /**  One of the list of Change Listeners */
   private ChangeListener cl;

   /**
    *   Fires a ChangeEvent from the WlzImgView.
    */
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
