/*!
* @file         WlzDBox2.java
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

public class WlzDBox2 extends WlzBase implements Cloneable
{
  // Version string.
  public static String ident = "Id$$";

  public double	xMin;
  public double	xMax;
  public double	yMin;
  public double	yMax;

  /*!
  * @brief:	Constructor
  * @param:     void
  */
  public	WlzDBox2()
  {
    xMin = 0.0;
    xMax = 0.0;
    yMin = 0.0;
    yMax = 0.0;
  }

  /*!
  * @brief    Constructor
  * @param     xmin		Minimum box x value
  * @param     ymin		Minimum box y value
  * @param     xmax		Maximum box x value
  * @param     ymax		Maximum box y value
  */
  public	WlzDBox2(double xmin, double ymin, double xmax, double ymax)
  {
    xMin = xmin;
    xMax = xmax;
    yMin = ymin;
    yMax = ymax;
  }

  /*!
  * @brief	Constructs a new Java (2D) rectangle which has the same
  *		coordinates as this Java Woolz box.
  */
  public Rectangle2D toRectangle()
  {
    return(new Rectangle2D.Double(xMin, yMin, xMax - xMin, yMax - yMin));
  }

  /*!
  * @return	Clone of this object.
  * \brief	Implements cloning
  */
  public Object	clone()
  {
    return(new WlzDBox2(xMin, yMin, xMax, yMax));
  }

  /*!
  * @return     true if this object is the same as the given object,
  *		otherwise false.
  * @brief    	Indicates whether some other object is "equal to" this
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
      			  ((WlzDBox2 )other).xMin) < Double.MIN_VALUE) &&
                (Math.abs(yMin -
			  ((WlzDBox2 )other).yMin) < Double.MIN_VALUE) &&
                (Math.abs(xMax -
			  ((WlzDBox2 )other).xMax) < Double.MIN_VALUE) &&
                (Math.abs(yMax -
			  ((WlzDBox2 )other).yMax) < Double.MIN_VALUE);
    }
    return(isEqual);
  }
}
