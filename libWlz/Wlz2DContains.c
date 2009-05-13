#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _Wlz2DContains_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/Wlz2DContains.c
* \author       Nick Burton
* \date         August 2002
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Takes a WLZ_2D_DOMAINOBJ, calls WlzLabel to split the domain
* 		and returns the one containing point(x,y).
* \ingroup	WlzBinaryOps
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <Wlz.h>

/*!
* \return	Object containing the point (x,y).
* \ingroup	WlzBinaryOps
* \brief	Takes a WLZ_2D_DOMAINOBJ, calls WlzLabel to split the domain
* 		and returns the one containing point(x,y).
* \param	obj		Given WLZ_2D_DOMAINOBJ object.
* \param	x		Column coordinate.
* \param	y		Line coordinate.
* \param	dstErr		Destination error code pointer, may be NULL.
*/
WlzObject *Wlz2DContains(WlzObject *obj,
                          double x,
                          double y,
                          WlzErrorNum *dstErr) {


   WlzObject    *retObj = NULL;
   WlzObject	**objArray;

   int          i, nobjs;
   /* int          maxNobjs=1024; */
   int          maxNobjs=2048;
   int 		found = 0;

   WlzErrorNum  errNum=WLZ_ERR_NONE;
   WlzConnectType connectivity = WLZ_8_CONNECTED;

   /*
   fprintf(stderr, "entering Wlz2DContains %f,%f\n", x,y);
   fflush(stderr);
   */

   if(obj->type != WLZ_2D_DOMAINOBJ) return (NULL);

   /* get array of domains */
   if(obj != NULL) {
      errNum = WlzLabel(obj,
                        &nobjs,
                        &objArray,
                        maxNobjs,
                        0,
                        connectivity);
   }
   /*
   fprintf(stderr, "got array of %d objects\n", nobjs);
   fflush(stderr);
   */

   /* select the required domain */
   /* ie the one which contains clicked point */

   if((errNum == WLZ_ERR_NONE) && (nobjs > 0)) {
      for(i=0; i<nobjs; i++) {
	 /*
	 fprintf(stderr, "checking domain # %d\n", i);
	 fflush(stderr);
	 */
	 if(WlzInsideDomain(objArray[i],
			    0.0,
			    y,
			    x,
			    &errNum)) {
	    if(errNum == WLZ_ERR_NONE) {
	       found = 1;
	       /*
	       fprintf(stderr, "domain # %d contains point\n", i);
	       fflush(stderr);
	       */
	       break;
	    } else {
	       (void )fprintf(stderr,
	                      "WlzInsideDomain, Wlz error: %d\n",
			      (int )errNum);
	       fflush(stderr);
	    }
	 } else {
	    WlzFreeObj(objArray[i]);
	 }
      }
      if((found == 1) && (errNum == WLZ_ERR_NONE)) {
	 /* retObj = objArray[i]; */
	 retObj = WlzCopyObject(objArray[i], &errNum);
	 WlzFreeObj(objArray[i]);
	 if(errNum != WLZ_ERR_NONE) {
	    (void )fprintf(stderr,
	                   "WlzCopyObject, Wlz error: %d\n",
			   (int )errNum);
	    fflush(stderr);
	 }
      }
   }

   *dstErr = errNum;
   /*
   fprintf(stderr, "leaving Wlz2DContains\n");
   fflush(stderr);
   */
   return retObj;
}
