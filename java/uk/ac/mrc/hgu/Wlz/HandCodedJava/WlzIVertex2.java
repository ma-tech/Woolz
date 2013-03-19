/*!
* @file         WlzIVertex2.java
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
* @brief        Java object to mirror the Woolz WlzIVertex2 structure.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;

import java.lang.*;
import java.awt.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzIVertex2 extends WlzBase implements Cloneable
{
  // Version string.
  public static String ident = "Id$$";

  public int	vtX;
  public int	vtY;

  /*!
  * @brief      Constructor
  */
  public	WlzIVertex2()
  {
    vtX = 0;
    vtY = 0;
  }

  /*!
  * @brief      Constructor
  * @param      x		The x coordinate.
  * @param      y		The y coordinate.
  */
  public	WlzIVertex2(int x, int y)
  {
    vtX = x;
    vtY = y;
  }

  /*!
  * @brief  	Constructs a new Java (2D) point which has the same
  *		coordinates as this Java Woolz vertex.
  */
  public Point toPoint()
  {
    return(new Point(vtX, vtY));
  }

  /*!
  * @return     Clone of this object.
  * @brief      Implements cloning.
  */
  public Object clone()
  {
    return(new WlzIVertex2(vtX, vtY));
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
      isEqual = (vtX == ((WlzIVertex2 )other).vtX) &&
                (vtY == ((WlzIVertex2 )other).vtY);
    }
    return(isEqual);
  }
}
