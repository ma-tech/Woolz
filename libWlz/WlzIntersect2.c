#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzIntersect2.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Tue Aug 19 08:08:06 2003
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzBinaryOps
* \brief        Calculates the intersection between two woolz
 domain objects
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <Wlz.h>

/* function:     WlzIntersect2    */
/*! 
* \ingroup      WlzBinaryOps
* \brief        Calculate the set intersection between two domain
 objects. This is a convenience routine calling WlzIntersectN with
 uvt=0. Input objects must be domain objects of the same type (2D
 or 3D) and non-NULL. Type WLZ_EMPTY_OBJ is legal, clearly an empty
 domain will be returned.
*
* \return       Intersection object with NULL value table
* \param    obj1	first input object
* \param    obj2	second input object
* \param    dstErr	error return.
* \par      Source:
*                WlzIntersect2.c
*/
WlzObject *WlzIntersect2(
  WlzObject *obj1,
  WlzObject *obj2,
  WlzErrorNum	*dstErr)
{
  WlzObject *objs[2];

  objs[0] = obj1;
  objs[1] = obj2;
  return(WlzIntersectN(2, objs, 0, dstErr));
}
