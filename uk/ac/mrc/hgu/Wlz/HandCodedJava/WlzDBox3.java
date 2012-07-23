/*!
* @file         WlzDBox3.java
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
* @brief        Java object to mirror the Woolz WlzDBox2 structure.
* @ingroup      JWlz
*/
package uk.ac.mrc.hgu.Wlz;

import java.lang.*;
import java.awt.geom.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzDBox3 extends WlzDBox2 implements Cloneable
{
  // Version string.
  public static String ident = "Id$$";

  public double	zMin;
  public double	zMax;

  /*!
  * @brief	Constructor
  */
  public	WlzDBox3()
  {
    zMin = 0.0;
    zMax = 0.0;
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
  public	WlzDBox3(double xmin, double ymin, double zmin,
  			 double xmax, double ymax, double zmax)
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
    return(new WlzDBox3(xMin, yMin, zMin, xMax, yMax, zMax));
  }

  /*!
  * @return     True if this object is the same as the given object,
  *		otherwise false.
  * @brief     	Indicates whether some other object is "equal to" this
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
      			  ((WlzDBox3 )other).xMin) < Double.MIN_VALUE) &&
                (Math.abs(yMin -
			  ((WlzDBox3 )other).yMin) < Double.MIN_VALUE) &&
                (Math.abs(zMin -
			  ((WlzDBox3 )other).zMin) < Double.MIN_VALUE) &&
                (Math.abs(xMax -
			  ((WlzDBox3 )other).xMax) < Double.MIN_VALUE) &&
                (Math.abs(yMax -
			  ((WlzDBox3 )other).yMax) < Double.MIN_VALUE) &&
                (Math.abs(zMax -
			  ((WlzDBox3 )other).zMax) < Double.MIN_VALUE);
    }
    return(isEqual);
  }
}
