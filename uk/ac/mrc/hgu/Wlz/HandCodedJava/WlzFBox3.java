/*!
* @file         WlzFBox3.java
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
* @brief        Java object to mirror the Woolz WlzFBox2 structure.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;

import java.lang.*;
import java.awt.geom.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzFBox3 extends WlzFBox2 implements Cloneable
{
  // Version string.
  public static String ident = "Id$$";

  public float	zMin;
  public float	zMax;

  /*!
  * @brief      Constructor
  * @param      void
  */
  public	WlzFBox3()
  {
    zMin = 0.0f;
    zMax = 0.0f;
  }

  /*!
  * @brief      Constructor
  * @param      xmin		Minimum box x value.
  * @param      ymin		Minimum box y value.
  * @param      zmin		Minimum box z value.
  * @param      xmax		Maximum box x value.
  * @param      ymax		Maximum box y value.
  * @param      zmax		Maximum box z value.
  */
  public	WlzFBox3(float xmin, float ymin, float zmin,
  			 float xmax, float ymax, float zmax)
  {
    xMin = xmin;
    xMax = xmax;
    yMin = ymin;
    yMax = ymax;
    zMin = zmin;
    zMax = zmax;
  }

  /*!
  * @return     Clone of this object.
  * @brief      Implements cloning.
  */
  public Object clone()
  {
    return(new WlzFBox3(xMin, yMin, zMin, xMax, yMax, zMax));
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
      isEqual = (Math.abs(xMin -
      			  ((WlzFBox3 )other).xMin) < Double.MIN_VALUE) &&
                (Math.abs(yMin -
			  ((WlzFBox3 )other).yMin) < Double.MIN_VALUE) &&
                (Math.abs(zMin -
			  ((WlzFBox3 )other).zMin) < Double.MIN_VALUE) &&
                (Math.abs(xMax -
			  ((WlzFBox3 )other).xMax) < Double.MIN_VALUE) &&
                (Math.abs(yMax -
			  ((WlzFBox3 )other).yMax) < Double.MIN_VALUE) &&
                (Math.abs(zMax -
			  ((WlzFBox3 )other).zMax) < Double.MIN_VALUE);
    }
    return(isEqual);
  }
}
