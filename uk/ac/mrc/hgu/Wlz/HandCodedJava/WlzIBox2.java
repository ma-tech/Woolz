/************************************************************************
* Project:      Java Woolz
* Title:        WlzIBox2.java
* Date:         January 1999
* Purpose:      Java object to mirror the Woolz WlzIBox2 structure.
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
import java.awt.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzIBox2 extends WlzBase implements Cloneable
{
  public int	xMin;
  public int	xMax;
  public int	yMin;
  public int	yMax;

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     void
  **********************************************************************/
  public	WlzIBox2()
  {
    xMin = 0;
    xMax = 0;
    yMin = 0;
    yMax = 0;
  }

  /**********************************************************************
  * Purpose:    Constructor
  * @param:     xmin		minimum box x value
  * @param:     ymin		minimum box y value
  * @param:     xmax		maximum box x value
  * @param:     ymax		maximum box y value
  **********************************************************************/
  public	WlzIBox2(int xmin, int ymin, int xmax, int ymax)
  {
    xMin = xmin;
    xMax = xmax;
    yMin = ymin;
    yMax = ymax;
  }

  /**********************************************************************
  * Purpose:	Constructs a new Java (2D) rectangle which has the same
  *		coordinates as this Java Woolz box.
  * @param:	void
  **********************************************************************/
  public Rectangle toRectangle()
  {
    return(new Rectangle(xMin, yMin, xMax - xMin + 1, yMax - yMin + 1));
  }

  /**********************************************************************
  * Purpose:    Implements cloning.
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzIBox2(xMin, yMin, xMax, yMax));
  }

  /**********************************************************************
  * Purpose:    Indicates whether some other object is "equal to" this
  *		Woolz pointer.
  * @return     true if this object is the same as the given object,
  *		otherwise false.
  * @param:     obj		the given object for comparison.
  **********************************************************************/
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
      isEqual = (xMin == ((WlzIBox2 )other).xMin) &&
                (yMin == ((WlzIBox2 )other).yMin) &&
                (xMax == ((WlzIBox2 )other).xMax) &&
                (yMax == ((WlzIBox2 )other).yMax);
    }
    return(isEqual);
  }
}
