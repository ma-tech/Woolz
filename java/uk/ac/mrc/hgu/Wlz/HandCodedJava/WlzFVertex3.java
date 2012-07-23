/*!
* @file         WlzFVertex3.java
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
* @brief        Java object to mirror the Woolz WlzFVertex3 structure.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;

import java.lang.*;
import java.awt.geom.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzFVertex3 extends WlzFVertex2 implements Cloneable
{
  // Version string.
  public static String ident = "Id$$";

  public float	vtZ;

  /*!
  * @brief      Constructor
  */
  public	WlzFVertex3()
  {
    vtZ = 0.0f;
  }

  /*!
  * @brief      Constructor
  * @param      x		The x coordinate.
  * @param      y		The y coordinate.
  * @param      z		The z coordinate.
  */
  public	WlzFVertex3(float x, float y, float z)
  {
    vtX = x;
    vtY = y;
    vtZ = z;
  }

  /*!
  * @return     Clone of this object.
  * @brief      Implements cloning.
  */
  public Object clone()
  {
    return(new WlzFVertex3(vtX, vtY, vtZ));
  }

  /*!
  * @return     true if this object is the same as the given object,
  *		otherwise false.
  * @brief      Indicates whether some other object is "equal to" this
  *		Woolz pointer.
  * @param      obj		The given object for comparison.
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
      isEqual = (Math.abs(vtX -
			  ((WlzFVertex3 )other).vtY) < Float.MIN_VALUE) &&
                (Math.abs(vtY -
			  ((WlzFVertex3 )other).vtY) < Float.MIN_VALUE) &&
                (Math.abs(vtZ -
			  ((WlzFVertex3 )other).vtZ) < Float.MIN_VALUE);
    }
    return(isEqual);
  }
}
