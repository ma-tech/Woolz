/************************************************************************
* Project:      Java Woolz
* Title:        WlzBase.java
* Date:         January 1999
* Purpose:      Base class for all Java Woolz classes.
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
import uk.ac.mrc.hgu.Wlz.*;

public class WlzBase extends Object
{
  /* Load the Woolz dynamic lib. */
  static
  {
    System.loadLibrary("JWlz");
  }
}
