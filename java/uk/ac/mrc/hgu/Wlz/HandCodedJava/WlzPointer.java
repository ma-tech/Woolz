/*!
* @file         WlzPointer.java
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
* @brief        Java binding for Woolz native pointers.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzPointer extends WlzNative
{
  // Version string.
  public static String ident = "Id$$";

  /*!
  * @return     True if this object is the same as the given object,
  *		otherwise false.
  * @brief      Indicates whether some other object is "equal to" this
  *		Woolz pointer.
  * @param      obj		the given object for comparison.
  */
  public boolean equals(Object other)
  {
    boolean	isEqual;

    if(other == null)
    {
      isEqual = false;
    }
    else if(other == this)
    {
      isEqual = true;
    }
    else if(WlzPointer.class != other.getClass())
    {
      isEqual = false;
    }
    else
    {
      isEqual = value == ((WlzPointer )other).value;
    }
    return(isEqual);
  }

  /*!
  * @return     Hash code value for the Woolz pointer.
  * @brief      Computes a hashcode for this pointer.
  */
  public int 	hashCode()
  {
    return((int )((value >>> 32) + (value & 0xFFFFFFFF)));
  }

  /*!
  * @return     True if the pointer is null.
  * @brief      Tests the pointer to see if it is null.
  */
  public boolean isNull()
  {
    return(value == 0);
  }

  /*!
  * @brief      Makes the pointer a null pointer.
  */
  public void	clear()
  {
    value = 0;
  }
}
