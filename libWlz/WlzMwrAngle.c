#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMwrangle.c
* \author       Jim Piper.
* \date         Fri Feb  8 13:12:02 2002
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2001 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzConvexHull
* \brief        Construct orientation of minimum width rectangle from convex hull of a woolz 2D domain-object.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <Wlz.h>

/*
 * extract angle (as scaled sin, cos) of minimum width rectangle
 * from convex hull
 * (it is relatively obvious that minimum width rectangle long
 * side must be parallel to a chord of convex hull, and all
 * sides must have at least one vertex of convex hull lying
 * within them).
 */


/* function:     WlzMwrAngle    */
/*! 
* \ingroup      WlzConvexHull
* \brief        extract angle (as scaled sin, cos) of minimum width rectangle
from convex hull (it is relatively obvious that minimum width rectangle long
 side must be parallel to a chord of convex hull, and all sides must have at
 least one vertex of convex hull lying within them).
*
* \return       angle of the minimum width rectangle
* \param    cvh	input convex hull
* \param    dstErr	destination error
* \par      Source:
*                WlzMwrangle.c
*/
double WlzMwrAngle(WlzObject *cvh, WlzErrorNum *dstErr)
{
  static double gap (WlzChord *ch, WlzPolygonDomain *pdom);
  WlzPolygonDomain *pdom;
  WlzChord *chr, *ch;
  WlzConvHullValues *cdom;
  int i;
  double c=1.0, s=0.0, h, w, minwidth;

  WlzErrorNum errNum = WLZ_ERR_NONE;

  if (cvh == NULL) {
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  else if (cvh->type == WLZ_EMPTY_OBJ) {
    return 0.0;
  }

  else if (cvh->type != WLZ_CONV_HULL) {
    errNum = WLZ_ERR_OBJECT_TYPE;
  }

  else {
    pdom = cvh->domain.poly;
    cdom = cvh->values.c;
    ch = cdom->ch;
    minwidth = gap(ch,pdom);
    chr = ch++;
    for (i=1; i<cdom->nchords; i++) {
      if ((w=gap(ch,pdom)) < minwidth) {
	minwidth = w;
	chr = ch;
      }
      ch++;
    }
    c = chr->acon;
    s = chr->bcon;
    h = sqrt (s*s + c*c);
    s /= h;
    c /= h;
  }

  if(dstErr)
  {
    *dstErr = errNum;
  }

  return(atan2(c,s));
}



/*
 * find width of polygon perpendicular to a chord.
 * also return the vertex numbers in polygon which are maximally distant.
 */
static double gap (WlzChord *ch, WlzPolygonDomain *pdom)
{
  WlzIVertex2 *vtx;
  int i;
  double c, minc, maxc, acon, bcon;
  vtx = pdom->vtx;
  acon = ch->acon;
  bcon = ch->bcon;
  for (i=0; i<pdom->nvertices; i++) {
    c = acon*vtx->vtX - bcon*vtx->vtY;
    if (i == 0) {
      minc = maxc = c;
    } else if (c > maxc) {
      maxc = c;
    } else if (c < minc) {
      minc = c;
    }
    vtx++;
  }
  /* return 32*gap (rounded to integer) */
  return((maxc - minc) / ch->cl);
}

