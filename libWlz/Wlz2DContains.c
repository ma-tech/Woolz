#pragma ident "MRC HGU $Id$"
/*!
* \file         Wlz2DContains.c
* \author       Nick Burton
* \date         August 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        takes a WLZ_2D_DOMAINOBJ, calls WlzLabel to split the domain
                and returns the one containing point(x,y)
* \todo         -
* \bug          None known.
*/
#include <stdio.h>
#include <Wlz.h>

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
   /* int hack = 1; */ /* for debugging with workshop */
   int hack = 0;

   while(hack)
   {
     sleep(1);
   }

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
   } // if
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
	       fprintf(stderr, "WlzInsideDomain, Wlz error: %d\n", errNum);
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
	    fprintf(stderr, "WlzCopyObject, Wlz error: %d\n", errNum);
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
