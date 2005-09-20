#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzDiffDomain.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Functions for computing the domain difference between
* 		objects.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

extern WlzObject *WlzDiffDomain3d(WlzObject *obj1,
				  WlzObject *obj2,
				  WlzErrorNum	*dstErr);

/*!
* \return	Object with domain equal to the set difference between the
*		first and second object, with valuetable from the first
*		object.
* \ingroup	WlzDomainOps
* \brief	Calculates the domain difference between two objects.
* \param	obj1			First object.
* \param	obj2			Second object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzDiffDomain(
  WlzObject *obj1,
  WlzObject *obj2,
  WlzErrorNum	*dstErr)
{
  WlzIntervalWSpace	iwsp2, iwsp1;
  WlzInterval		*intp;
  WlzInterval 		*jntp;
  WlzDomain 		diffdom;
  WlzValues		values;
  WlzIntervalDomain 	*idom1;
  int 			k1, dfinished, intervalcount;
  WlzObject 		*diff=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object 1 */
  if( obj1 == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj1->type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      if( obj1->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_EMPTY_OBJ:
      diffdom.i = NULL;
      values.v = NULL;
      return WlzMakeMain(WLZ_EMPTY_OBJ, diffdom, values, NULL,
			 NULL, dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    }
  }

  /* check object 2 - obj2 == NULL is now an error */
  if( (errNum == WLZ_ERR_NONE) && (obj2 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj2->type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      if( obj2->type != obj1->type ){
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
      if( obj2->domain.core == NULL ){
	errNum = WLZ_ERR_OBJECT_NULL;
	break;
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeMain(obj1->type, obj1->domain, obj1->values,
			 NULL, NULL, dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;

    }
  }

  /* switch for 3D objects */
  if( (errNum == WLZ_ERR_NONE) && (obj1->type == WLZ_3D_DOMAINOBJ) ){
    return WlzDiffDomain3d( obj1, obj2, dstErr );
  }

  /*
   * make a new interval table that is big enough
   */
  if( errNum == WLZ_ERR_NONE ){
    idom1 = obj1->domain.i;
    k1 = idom1->kol1;
    if( (intp = (WlzInterval *)
	 AlcCalloc(WlzIntervalCount(idom1, NULL)+
		   WlzIntervalCount(obj2->domain.i, NULL),
		   sizeof(WlzInterval))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    jntp = intp;
  }

  if( errNum == WLZ_ERR_NONE ){
    if( (diffdom.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					   idom1->line1, idom1->lastln,
					   k1, idom1->lastkl,
					   &errNum)) == NULL ){
      AlcFree((void *) intp);
    }
    else {
      diffdom.i->freeptr = AlcFreeStackPush(diffdom.i->freeptr, (void *)intp,
					    NULL);
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    if( (diff = WlzMakeMain(WLZ_2D_DOMAINOBJ, diffdom, obj1->values,
			    NULL, NULL, &errNum)) == NULL ){
      WlzFreeIntervalDomain(diffdom.i);
    }
  }

  /*
   * initialise interval scanning
   */
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitRasterScan(obj1, &iwsp1, WLZ_RASTERDIR_ILIC);
  }
  if( errNum == WLZ_ERR_NONE ){
    errNum = WlzInitRasterScan(obj2, &iwsp2, WLZ_RASTERDIR_ILIC);
    dfinished = 0;
    intervalcount = 0;
  }

  /*
   * scan object constructing intervals after subtraction of obj2
   */
  if( (errNum == WLZ_ERR_NONE) &&
     ((errNum = WlzNextInterval(&iwsp2)) == WLZ_ERR_NONE) ){

    while( (errNum = WlzNextInterval(&iwsp1)) == WLZ_ERR_NONE ){
      /*
       * case 1 - interval in obj1 succedes interval in obj2 -
       *	    get next interval of obj2
       */
      while( (errNum == WLZ_ERR_NONE) && (dfinished == 0) &&
	     (iwsp1.linpos > iwsp2.linpos ||
	      (iwsp1.linpos == iwsp2.linpos &&
	       iwsp1.lftpos > iwsp2.rgtpos))){
	switch( errNum = WlzNextInterval(&iwsp2) ){
	case WLZ_ERR_NONE:
	  dfinished = 0;
	  break;
	case WLZ_ERR_EOO:
	  errNum = WLZ_ERR_NONE;
	  dfinished = 1;
	  break;
	default:
	  break;
	}
      }
      if( errNum != WLZ_ERR_NONE ){
	break;
      }
      /*
     * case 2 - interval in obj1 precedes interval in obj2
     *	    or obj2 finished - just copy interval
     */
      if( (dfinished != 0) ||
	 (iwsp1.linpos < iwsp2.linpos ||
	  ((iwsp1.linpos == iwsp2.linpos) &&
	   (iwsp1.rgtpos < iwsp2.lftpos))
	   )
	){
	jntp->ileft = iwsp1.lftpos - k1;
	jntp->iright = iwsp1.rgtpos - k1;
	intervalcount++;
	jntp++;
      }
      /*
     * case 3 - intervals overlap.
     * there may be more than one obj2 interval.
     */
      else {
	/*
	 * test left end of obj1 interval < left
	 * end of obj2 interval, when complete
	 * interval can be constructed immediately
	 */
	if (iwsp2.lftpos > iwsp1.lftpos) {
	  jntp->ileft = iwsp1.lftpos - k1;
	  jntp->iright = iwsp2.lftpos -1 - k1;
	  intervalcount++;
	  jntp++;
	}
	do {
	  /*
	   * if right end of obj2 interval <
	   * right end of obj1 interval can
	   * plant left end of new interval
	   */
	  if (iwsp2.rgtpos < iwsp1.rgtpos) {
	    jntp->ileft = iwsp2.rgtpos +1 - k1;
	    switch( errNum = WlzNextInterval(&iwsp2) ){
	    case WLZ_ERR_NONE:
	      dfinished = 0;
	      break;
	    case WLZ_ERR_EOO:
	      errNum = WLZ_ERR_NONE;
	      dfinished = 1;
	      break;
	    default:
	      break;
	    }
	    if( errNum != WLZ_ERR_NONE ){
	      break;
	    }
	    
	    /*
	     * check if next obj2 interval
	     * still overlaps this obj1 interval
	     */
	    if( (dfinished == 0) &&
	       (iwsp1.linpos == iwsp2.linpos) &&
	       (iwsp1.rgtpos >= iwsp2.lftpos))
	      jntp->iright = iwsp2.lftpos -1 -k1;
	    else
	      jntp->iright = iwsp1.rgtpos - k1;
	    intervalcount++;
	    jntp++;
	  }
	} while( (dfinished == 0) && 
		(iwsp1.linpos == iwsp2.linpos) &&
		(iwsp1.rgtpos > iwsp2.rgtpos) );
      }
      if( errNum != WLZ_ERR_NONE ){
	break;
      }

      /*
     * is this obj1 interval the last in the line ?
     * if so, load up the new intervals into the
     * diff domain.
     */
      if (iwsp1.intrmn == 0) {
	errNum = WlzMakeInterval(iwsp1.linpos, diffdom.i,
				 intervalcount, intp);
	intervalcount = 0;
	intp = jntp;
      }
      if( errNum != WLZ_ERR_NONE ){
	break;
      }
    }
  }
  if(errNum == WLZ_ERR_EOO) 	        /* Reset error from end of intervals */ 
  {
    errNum = WLZ_ERR_NONE;
  }

  /* this checks the area and returns NULL if zero
     in principle it should return an "empty object" */
  if( (errNum == WLZ_ERR_NONE) && !WlzIntervalCount(diffdom.i, NULL) ){
    WlzFreeObj( diff );
    diffdom.i = NULL;
    values.v = NULL;
    diff = WlzMakeMain(WLZ_EMPTY_OBJ, diffdom, values,
			NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return diff;
}
