#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzMwrAngle_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzMwrAngle.c
* \author       Jim Piper
* \date         February 2002
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
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
* \brief	Computes the minimum width rectangle from a
* 		convex hull.
* \ingroup	WlzConvexHull
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

static double gap (WlzChord *ch, WlzPolygonDomain *pdom);

/* function:     WlzMwrAngle    */
/*! 
* \ingroup      WlzConvexHull
* \return       angle of the minimum width rectangle
* \brief        extract angle (as scaled sin, cos) of minimum width rectangle
*		from convex hull (it is relatively obvious that minimum width
*		rectangle long side must be parallel to a chord of convex hull,
*		and all sides must have at least one vertex of convex hull
*		lying within them).
*
* \param    cvh	input convex hull
* \param    dstErr	destination error
*/
double WlzMwrAngle(WlzObject *cvh, WlzErrorNum *dstErr)
{
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
