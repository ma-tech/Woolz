/************************************************************************
* Project:      Java Woolz
* Title:        WlzPixelV.java
* Date:         January 1999
* Purpose:      Java binding for Woolz grey value structure.
* Copyright:	1997 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Maintenance:	Log changes below, with most recent at top of list.
* @author       Bill Hill (bill@hgu.mrc.ac.uk)
* @version 	MRC HGU %I%, %G%
************************************************************************/
package uk.ac.mrc.hgu.Wlz;
import java.lang.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzPixelV extends WlzNative implements Cloneable
{
  protected int	type;

  /**********************************************************************
  * Purpose:    Constructor for use only by C side of the JNI.
  * @param:     t			the type of pixel
  *		v			the pixel value.
  **********************************************************************/
  private WlzPixelV(int t, long v)
  {
    type = t;
    value = v;
  }

  /**********************************************************************
  * Purpose:    Constructor for a long valued pixel.
  * @param:     v			the pixel value.
  **********************************************************************/
  public	WlzPixelV(long v)
  {
    type = WlzGreyType.WLZ_GREY_LONG;
    setLongValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long valued pixel.
  * @param:     v			the pixel value.
  **********************************************************************/
  public	WlzPixelV(int v)
  {
    type = WlzGreyType.WLZ_GREY_INT;
    setIntValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long valued pixel.
  * @param:     v			the pixel value.
  **********************************************************************/
  public	WlzPixelV(short v)
  {
    type = WlzGreyType.WLZ_GREY_SHORT;
    setShortValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long valued pixel.
  * @param:     v			the pixel value.
  **********************************************************************/
  public	WlzPixelV(byte v)
  {
    type = WlzGreyType.WLZ_GREY_UBYTE;
    setByteValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long valued pixel.
  * @param:     v			the pixel value.
  **********************************************************************/
  public	WlzPixelV(float v)
  {
    type = WlzGreyType.WLZ_GREY_FLOAT;
    setFloatValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long valued pixel.
  * @param:     v			the pixel value.
  **********************************************************************/
  public	WlzPixelV(double v)
  {
    type = WlzGreyType.WLZ_GREY_DOUBLE;
    setDoubleValue(v);
  }

  /**********************************************************************
  * Purpose:    Get the grey value type.
  * @return     int		the grey value type of the pixel
  * @param:     void
  **********************************************************************/
  public int	getType()
  {
    return((int )type);
  }

  /**********************************************************************
  * Purpose:    Get the long grey value.
  * @return     void
  * @param:     void
  **********************************************************************/
  public long	getLongValue()
  {
    long	val = 0;

    switch(type)
    {
      case WlzGreyType.WLZ_GREY_LONG:
        val = super.getLongValue();
	break;
      case WlzGreyType.WLZ_GREY_INT:
	val = (long )(super.getIntValue());
	break;
      case WlzGreyType.WLZ_GREY_SHORT:
	val = (long )(super.getShortValue());
	break;
      case WlzGreyType.WLZ_GREY_UBYTE:
	val = (long )(super.getByteValue());
	break;
      case WlzGreyType.WLZ_GREY_FLOAT:
	val = (long )(Math.round(super.getFloatValue()));
	break;
      case WlzGreyType.WLZ_GREY_DOUBLE:
	val = Math.round(super.getDoubleValue());
        break;
    }
    return(val);
  }

  /**********************************************************************
  * Purpose:    Get the int grey value.
  * @return     void
  * @param:     void
  **********************************************************************/
  public int	getIntValue()
  {
    int		val = 0;

    switch(type)
    {
      case WlzGreyType.WLZ_GREY_LONG:
        val = (int )(super.getLongValue());
	break;
      case WlzGreyType.WLZ_GREY_INT:
	val = super.getIntValue();
	break;
      case WlzGreyType.WLZ_GREY_SHORT:
	val = (int )(super.getShortValue());
	break;
      case WlzGreyType.WLZ_GREY_UBYTE:
	val = (int )(super.getByteValue());
	break;
      case WlzGreyType.WLZ_GREY_FLOAT:
	val = Math.round(super.getFloatValue());
	break;
      case WlzGreyType.WLZ_GREY_DOUBLE:
	val = (int )(Math.round(super.getDoubleValue()));
        break;
    }
    return(val);
  }

  /**********************************************************************
  * Purpose:    Get the short grey value.
  * @return     void
  * @param:     void
  **********************************************************************/
  public short	getShortValue()
  {
    short	val = 0;

    switch(type)
    {
      case WlzGreyType.WLZ_GREY_LONG:
        val = (short )(super.getLongValue());
	break;
      case WlzGreyType.WLZ_GREY_INT:
	val = (short )(super.getIntValue());
	break;
      case WlzGreyType.WLZ_GREY_SHORT:
	val = super.getShortValue();
	break;
      case WlzGreyType.WLZ_GREY_UBYTE:
	val = (short )(super.getByteValue());
	break;
      case WlzGreyType.WLZ_GREY_FLOAT:
	val = (short )(Math.round(super.getFloatValue()));
	break;
      case WlzGreyType.WLZ_GREY_DOUBLE:
	val = (short )(Math.round(super.getDoubleValue()));
        break;
    }
    return(val);
  }

  /**********************************************************************
  * Purpose:    Get the byte grey value.
  * @return     void
  * @param:     void
  **********************************************************************/
  public byte	getByteValue()
  {
    byte	val = 0;

    switch(type)
    {
      case WlzGreyType.WLZ_GREY_LONG:
        val = (byte )(super.getLongValue());
	break;
      case WlzGreyType.WLZ_GREY_INT:
	val = (byte )(super.getIntValue());
	break;
      case WlzGreyType.WLZ_GREY_SHORT:
	val = (byte )(super.getShortValue());
	break;
      case WlzGreyType.WLZ_GREY_UBYTE:
	val = super.getByteValue();
	break;
      case WlzGreyType.WLZ_GREY_FLOAT:
	val = (byte )(Math.round(super.getFloatValue()));
	break;
      case WlzGreyType.WLZ_GREY_DOUBLE:
	val = (byte )(Math.round(super.getDoubleValue()));
        break;
    }
    return(val);
  }

  /**********************************************************************
  * Purpose:    Get the float grey value.
  * @return     void
  * @param:     void
  **********************************************************************/
  public float	getFloatValue()
  {
    float	val = 0.0f;

    switch(type)
    {
      case WlzGreyType.WLZ_GREY_LONG:
        val = (float )(super.getLongValue());
	break;
      case WlzGreyType.WLZ_GREY_INT:
	val = (float )(super.getIntValue());
	break;
      case WlzGreyType.WLZ_GREY_SHORT:
	val = (float )(super.getShortValue());
	break;
      case WlzGreyType.WLZ_GREY_UBYTE:
	val = (float )(super.getByteValue());
	break;
      case WlzGreyType.WLZ_GREY_FLOAT:
	val = super.getFloatValue();
	break;
      case WlzGreyType.WLZ_GREY_DOUBLE:
	val = (float )(super.getDoubleValue());
        break;
    }
    return(val);
  }

  /**********************************************************************
  * Purpose:    Get the double grey value.
  * @return     void
  * @param:     void
  **********************************************************************/
  public double	getDoubleValue()
  {
    double	val = 0.0;

    switch(type)
    {
      case WlzGreyType.WLZ_GREY_LONG:
        val = (double )(super.getLongValue());
	break;
      case WlzGreyType.WLZ_GREY_INT:
	val = (double )(super.getIntValue());
	break;
      case WlzGreyType.WLZ_GREY_SHORT:
	val = (double )(super.getShortValue());
	break;
      case WlzGreyType.WLZ_GREY_UBYTE:
	val = (double )(super.getByteValue());
	break;
      case WlzGreyType.WLZ_GREY_FLOAT:
	val = (double )(super.getFloatValue());
	break;
      case WlzGreyType.WLZ_GREY_DOUBLE:
	val = super.getDoubleValue();
        break;
    }
    return(val);
  }

  /**********************************************************************
  * Purpose:    Set the long grey value.
  * @return     void
  * @param:     val		the long value
  **********************************************************************/
  public void	setValue(long val)
  {
    type = WlzGreyType.WLZ_GREY_LONG;
    setLongValue(val);
  }

  /**********************************************************************
  * Purpose:    Set the int grey value.
  * @return     void
  * @param:     val		the int value
  **********************************************************************/
  public void	setValue(int val)
  {
    type = WlzGreyType.WLZ_GREY_INT;
    setIntValue(val);
  }

  /**********************************************************************
  * Purpose:    Set the short grey value.
  * @return     void
  * @param:     val		the short value
  **********************************************************************/
  public void	setValue(short val)
  {
    type = WlzGreyType.WLZ_GREY_SHORT;
    setShortValue(val);
  }

  /**********************************************************************
  * Purpose:    Set the byte grey value.
  * @return     void
  * @param:     val		the byte value
  **********************************************************************/
  public void	setValue(byte val)
  {
    type = WlzGreyType.WLZ_GREY_UBYTE;
    setByteValue(val);
  }

  /**********************************************************************
  * Purpose:    Set the float grey value.
  * @return     void
  * @param:     val		the float value
  **********************************************************************/
  public void	setValue(float val)
  {
    type = WlzGreyType.WLZ_GREY_FLOAT;
    setFloatValue(val);
  }

  /**********************************************************************
  * Purpose:    Set the double grey value.
  * @return     void
  * @param:     val		the double value
  **********************************************************************/
  public void	setValue(double val)
  {
    type = WlzGreyType.WLZ_GREY_DOUBLE;
    setDoubleValue(val);
  }

  /**********************************************************************
  * Purpose:    Implements cloning.
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzPixelV(type, value));
  }
}
