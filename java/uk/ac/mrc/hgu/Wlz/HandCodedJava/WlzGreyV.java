/************************************************************************
* Project:      Java Woolz
* Title:        WlzGreyV.java
* Date:         January 1999
* Purpose:      Java binding for Woolz grey value structure.
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
import uk.ac.mrc.hgu.Wlz.*;

public class WlzGreyV extends WlzNative implements Cloneable
{
  /**********************************************************************
  * Purpose:    Constructor for a long grey value.
  * @param:     v			the grey value.
  **********************************************************************/
  public	WlzGreyV(long v)
  {
    setLongValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long grey value.
  * @param:     v			the grey value.
  **********************************************************************/
  public	WlzGreyV(int v)
  {
    setIntValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long grey value.
  * @param:     v			the grey value.
  **********************************************************************/
  public	WlzGreyV(short v)
  {
    setShortValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long grey value.
  * @param:     v			the grey value.
  **********************************************************************/
  public	WlzGreyV(byte v)
  {
    setByteValue(v);
  }

  /**********************************************************************
  * Purpose:    Constructor for a long grey value.
  * @param:     v			the grey value.
  **********************************************************************/
  public	WlzGreyV(float v)
  {
    setFloatValue(v);
  }

  /**********************************************************************
  * Purpose:    Implements cloning.
  * @return:    Clone of this object.
  * @param:     void
  **********************************************************************/
  public Object clone()
  {
    return(new WlzGreyV(value));
  }

}
