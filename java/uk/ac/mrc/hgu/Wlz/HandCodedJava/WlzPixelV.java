/*!
* @file         WlzPixelV.java
* @author       Bill Hill
* @date         January 1999
* @version      $Id$
* @par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* @par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* @brief        Java binding for Woolz grey value structure.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;
import java.lang.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzPixelV extends WlzNative implements Cloneable
{
  // Version string.
  public static String ident = "Id$$";

  protected int	type;

  /*!
  * @brief      Constructor for use only by C side of the JNI.
  * @param      t			The type of pixel.
  * @param	v			The pixel value.
  */
  private WlzPixelV(int t, long v)
  {
    type = t;
    value = v;
  }

  /*!
  * @brief      Constructor for a long valued pixel.
  * @param      v			The pixel value.
  */
  public	WlzPixelV(long v)
  {
    type = WlzGreyType.WLZ_GREY_LONG;
    setLongValue(v);
  }

  /*!
  * @brief      Constructor for a long valued pixel.
  * @param      v			The pixel value.
  */
  public	WlzPixelV(int v)
  {
    type = WlzGreyType.WLZ_GREY_INT;
    setIntValue(v);
  }

  /*!
  * @brief      Constructor for a long valued pixel.
  * @param      v			The pixel value.
  */
  public	WlzPixelV(short v)
  {
    type = WlzGreyType.WLZ_GREY_SHORT;
    setShortValue(v);
  }

  /*!
  * @brief      Constructor for a long valued pixel.
  * @param      v			The pixel value.
  */
  public	WlzPixelV(byte v)
  {
    type = WlzGreyType.WLZ_GREY_UBYTE;
    setByteValue(v);
  }

  /*!
  * @brief      Constructor for a long valued pixel.
  * @param      v			The pixel value.
  */
  public	WlzPixelV(float v)
  {
    type = WlzGreyType.WLZ_GREY_FLOAT;
    setFloatValue(v);
  }

  /*!
  * @brief      Constructor for a long valued pixel.
  * @param      v			The pixel value.
  */
  public	WlzPixelV(double v)
  {
    type = WlzGreyType.WLZ_GREY_DOUBLE;
    setDoubleValue(v);
  }

  /*!
  * @return     int		The grey value type of the pixel.
  * @brief      Get the grey value type.
  */
  public int	getType()
  {
    return((int )type);
  }

  /*!
  * @brief      Get the long grey value.
  */
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

  /*!
  * @brief      Get the int grey value.
  */
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

  /*!
  * @brief      Get the short grey value.
  */
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

  /*!
  * @brief      Get the byte grey value.
  */
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

  /*!
  * @brief      Get the float grey value.
  */
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

  /*!
  * @brief      Get the double grey value.
  */
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

  /*!
  * @brief      Set the long grey value.
  * @param      val		the long value
  */
  public void	setValue(long val)
  {
    type = WlzGreyType.WLZ_GREY_LONG;
    setLongValue(val);
  }

  /*!
  * @brief      Set the int grey value.
  * @param      val		the int value
  */
  public void	setValue(int val)
  {
    type = WlzGreyType.WLZ_GREY_INT;
    setIntValue(val);
  }

  /*!
  * @brief      Set the short grey value.
  * @param      val		the short value
  */
  public void	setValue(short val)
  {
    type = WlzGreyType.WLZ_GREY_SHORT;
    setShortValue(val);
  }

  /*!
  * @brief      Set the byte grey value.
  * @param      val		the byte value
  */
  public void	setValue(byte val)
  {
    type = WlzGreyType.WLZ_GREY_UBYTE;
    setByteValue(val);
  }

  /*!
  * @brief      Set the float grey value.
  * @param      val		the float value
  */
  public void	setValue(float val)
  {
    type = WlzGreyType.WLZ_GREY_FLOAT;
    setFloatValue(val);
  }

  /*!
  * @brief      Set the double grey value.
  * @param      val		the double value
  */
  public void	setValue(double val)
  {
    type = WlzGreyType.WLZ_GREY_DOUBLE;
    setDoubleValue(val);
  }

  /*!
  * @brief      Implements cloning.
  * @return     Clone of this object.
  */
  public Object clone()
  {
    return(new WlzPixelV(type, value));
  }
}
