#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzUnion2.c
* \author       richard <Richard.Baldock@hgu.mrc.ac.uk>
* \date         Tue Aug 19 08:26:52 2003
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
* \brief        Convenience procedure to calculate the union of two
 woolz domain objects.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdio.h>
#include <Wlz.h>

/* function:     WlzUnion2    */
/*! 
* \ingroup      WlzBinaryOps
* \brief        Convenience procedure to calculate the union of two
 woolz domain objects. This calls WlzUnnionN() with uvt=0. Objects must
 be of the same type.
*
* \return       Union object pointer.
* \param    obj1	First inout object.
* \param    obj2	Second input object.
* \param    dstErr	errro return.
* \par      Source:
*                WlzUnion2.c
*/
WlzObject *WlzUnion2(WlzObject *obj1,
		     WlzObject *obj2,
		     WlzErrorNum *dstErr)
{
  WlzObject *objs[2];

  objs[0] = obj1;
  objs[1] = obj2;

  return WlzUnionN(2, objs, 0, dstErr);
}
