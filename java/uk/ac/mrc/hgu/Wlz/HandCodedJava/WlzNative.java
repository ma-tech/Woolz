/*!
* @file         WlzNative.java
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
* @brief        Java binding for Woolz native types; ie: int, short, byte,
* 		float, double and pointer.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzNative extends WlzBase
{
  // Version string.
  public static String ident = "Id$$";

  // Possible types held in the native value union.
  public static final int LONG = 0;
  public static final int INT = 1;
  public static final int SHORT = 2;
  public static final int BYTE = 3;
  public static final int FLOAT = 4;
  public static final int DOUBLE = 5;
  public static final int POINTER = 6;

  // The native value union.
  protected long value;

  // Native methods.
  private native long JWlzNativeUnionFromLong(long val);
  private native long JWlzNativeUnionFromInt(int val);
  private native long JWlzNativeUnionFromShort(short val);
  private native long JWlzNativeUnionFromByte(byte val);
  private native long JWlzNativeUnionFromFloat(float val);
  private native long JWlzNativeUnionFromDouble(double val);
  private native long JWlzNativeUnionToLong(long value);
  private native int JWlzNativeUnionToInt(long value);
  private native short JWlzNativeUnionToShort(long value);
  private native byte JWlzNativeUnionToByte(long value);
  private native float JWlzNativeUnionToFloat(long value);
  private native double JWlzNativeUnionToDouble(long value);

  /*!
  * @brief      Constructor, just to make it explicit.
  */
  public WlzNative()
  {
    value = 0;
  }

  /*!
  * @brief      Constructor with initial value which is only seen from
  *		the JNI C side.
  * @param      val			Initial value.
  */
  private WlzNative(long v)
  {
    value = v;
  }

  /*!
  * @brief      Set the value to the given long value.
  * @param      val			Given value.
  */
  protected void setLongValue(long val)
  {
    value = val;
    // value = JWlzNativeUnionFromLong(val);
  }

  /*!
  * @brief      Set the value to the given int value.
  * @param      val			Given value.
  */
  protected void setIntValue(int val)
  {
    value = JWlzNativeUnionFromInt(val);
  }

  /*!
  * @brief      Set the value to the given short value.
  * @param      val			Given value.
  */
  protected void setShortValue(short val)
  {
    value = JWlzNativeUnionFromShort(val);
  }

  /*!
  * @brief      Set the value to the given byte value.
  * @param      val			Given value.
  */
  protected void setByteValue(byte val)
  {
    value = JWlzNativeUnionFromByte(val);
  }

  /*!
  * @brief      Set the value to the given float value.
  * @param      val			Given value.
  */
  protected void setFloatValue(float val)
  {
    value = JWlzNativeUnionFromFloat(val);
  }

  /*!
  * @brief      Set the value to the given double value.
  * @param      val			Given value.
  */
  protected void setDoubleValue(double val)
  {
    value = JWlzNativeUnionFromDouble(val);
  }

  /*!
  * @brief      Get the long value.
  */
  protected long getLongValue()
  {
    // return(JWlzNativeUnionToLong(value));
    return(value);
  }

  /*!
  * @brief      Get the int value.
  */
  protected int getIntValue()
  {
    return(JWlzNativeUnionToInt(value));
  }

  /*!
  * @brief      Get the short value.
  */
  protected short getShortValue()
  {
    return(JWlzNativeUnionToShort(value));
  }

  /*!
  * @brief      Get the byte value.
  */
  protected byte getByteValue()
  {
    return(JWlzNativeUnionToByte(value));
  }

  /*!
  * @brief      Get the float value.
  */
  protected float getFloatValue()
  {
    return(JWlzNativeUnionToFloat(value));
  }

  /*!
  * @brief      Get the double value.
  */
  protected double getDoubleValue()
  {
    return(JWlzNativeUnionToDouble(value));
  }
}
