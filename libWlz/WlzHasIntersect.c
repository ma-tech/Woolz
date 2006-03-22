#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzHasIntersect_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzHasIntersect.c
* \author       Richard Baldock
* \date         December 2000
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
* \brief	Determines whether two objects have an intersection.
* \ingroup	WlzDomainUtils
* \todo         -
* \bug          None known.
*/


#include <stdlib.h>
#include <Wlz.h>


/*!
* \ingroup     Wlz
*/

static int WlzIntervalHasIntersect(
  WlzInterval	*intvl1,
  int		off1,
  WlzInterval	*intvl2,
  int		off2)
{
  if( (intvl1->iright + off1) < (intvl2->ileft + off2) ){
    return -1;
  }
  else if( (intvl2->iright + off2) < (intvl1->ileft + off1) ){
    return 1;
  }
  else {
    return 0;
  }
}

static int WlzIntervalLineHasIntersect(
  WlzIntervalLine	*intvln1,
  int			off1,
  WlzIntervalLine	*intvln2,
  int			off2)
{
  WlzInterval	*intvls1=intvln1->intvs;
  WlzInterval	*intvls2=intvln2->intvs;
  int	nintvs1=intvln1->nintvs;
  int	nintvs2=intvln2->nintvs;
  int	finished=0;

  while( !finished && nintvs1 && nintvs2 ){
    /* check current intervals */
    switch( WlzIntervalHasIntersect(intvls1, off1, intvls2, off2) ){
    case -1:
      intvls1++;
      nintvs1--;
      break;

    case 0:
      return 1;

    case 1:
      intvls2++;
      nintvs2--;
      break;
    }
  }
  return 0;
}


/* function:     WlzHasIntersection    */
/*! 
* \ingroup      WlzDomainUtils
* \brief        Determine if two domain objects intersect. The objects
must WLZ_2D_DOMAINOBJ, WLZ_3D_DOMAINOBJ or WLZ_EMPTY_OBJ. If neither object
is empty then they must be of the same type.
*
* \return       non-zero if the input objects intersect
* \param    obj1	input object 1
* \param    obj2	input object 2
* \param    dstErr	error return
* \par      Source:
*                WlzHasIntersect.c
*/
int WlzHasIntersection(
  WlzObject	*obj1,
  WlzObject	*obj2,
  WlzErrorNum	*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  int		rtnVal=0;
  int		i, p, pMin, pMax;
  WlzObject	*tmpObj1, *tmpObj2;
  WlzValues	values;

  /* check objects */
  if( (obj1 == NULL) || (obj2 == NULL) ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj1->type ){
    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      if( obj1->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_EMPTY_OBJ:
      break;

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }

    if( errNum == WLZ_ERR_NONE ){
      switch( obj2->type ){
      case WLZ_2D_DOMAINOBJ:
      case WLZ_3D_DOMAINOBJ:
	if( obj2->domain.core == NULL ){
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	break;

      case WLZ_EMPTY_OBJ:
	break;

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }

    if( errNum == WLZ_ERR_NONE ){
      if( obj1->type != obj2->type ){
	errNum = WLZ_ERR_OBJECT_TYPE;
      }
      else if(WlzIsEmpty(obj1, &errNum) ||
	      WlzIsEmpty(obj2, &errNum)){
	if( dstErr ){
	  *dstErr = errNum;
	}
	return 0;
      }
    }
  }

  /* valid objects now check for intersection */
  if( errNum == WLZ_ERR_NONE ){
    switch( obj1->type ){
    case WLZ_2D_DOMAINOBJ:
      /* check bounding box */
      if((obj1->domain.i->kol1 > obj2->domain.i->lastkl) ||
	 (obj2->domain.i->kol1 > obj1->domain.i->lastkl) ||
	 (obj1->domain.i->line1 > obj2->domain.i->lastln) ||
	 (obj2->domain.i->line1 > obj1->domain.i->lastln)){
	if( dstErr ){
	  *dstErr = errNum;
	}
	return 0;
      }

      /* if either rectangular then there is intersection */
      if((obj1->domain.i->type == WLZ_INTERVALDOMAIN_RECT) ||
	 (obj2->domain.i->type == WLZ_INTERVALDOMAIN_RECT)){
	rtnVal = 1;
      }
      /* iterate over common lines checking intervals */
      else {
	WlzIntervalWSpace	iwsp1, iwsp2;
	int			lMin, lMax;

	if( (errNum = WlzInitRasterScan(obj1, &iwsp1, WLZ_RASTERDIR_ILIC))
	   == WLZ_ERR_NONE ){
	  errNum = WlzInitRasterScan(obj2, &iwsp2, WLZ_RASTERDIR_ILIC);
	}

	lMin = WLZ_MAX(obj1->domain.i->line1, obj2->domain.i->line1);
	lMax = WLZ_MIN(obj1->domain.i->lastln, obj2->domain.i->lastln);

	for(i=lMin; i <= lMax; i++){
	  while((errNum == WLZ_ERR_NONE) &&
		(iwsp1.linpos < i)){
	    WlzNextInterval(&iwsp1);
	  }
	  if( errNum == WLZ_ERR_EOO ){
	    errNum = WLZ_ERR_NONE;
	    break;
	  }

	  while((errNum == WLZ_ERR_NONE) &&
		(iwsp2.linpos < i)){
	    WlzNextInterval(&iwsp2);
	  }
	  if( errNum == WLZ_ERR_EOO ){
	    errNum = WLZ_ERR_NONE;
	    break;
	  }

	  if( (iwsp1.linpos != i) || (iwsp2.linpos != i) ){
	    continue;
	  }
	
	  if( WlzIntervalLineHasIntersect(iwsp1.intvln, iwsp1.intdmn->kol1,
					  iwsp2.intvln, iwsp2.intdmn->kol1) ){
	    rtnVal = 1;
	    break;
	  }
	}
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      /* get plane bounds */
      pMin = WLZ_MAX(obj1->domain.p->plane1, obj2->domain.p->plane1);
      pMax = WLZ_MIN(obj1->domain.p->lastpl, obj2->domain.p->lastpl);
      values.core = NULL;
      for(p=pMin;(errNum == WLZ_ERR_NONE)&&(rtnVal == 0)&&(p <= pMax); p++){
	if( obj1->domain.p->domains[p - obj1->domain.p->plane1].core == NULL ){
	  continue;
	}
	if( obj2->domain.p->domains[p - obj2->domain.p->plane1].core == NULL ){
	  continue;
	}
	tmpObj1 = 
	  WlzMakeMain(WLZ_2D_DOMAINOBJ, 
		      obj1->domain.p->domains[p - obj1->domain.p->plane1],
		      values, NULL, NULL, NULL);
	tmpObj2 = 
	  WlzMakeMain(WLZ_2D_DOMAINOBJ, 
		      obj2->domain.p->domains[p - obj2->domain.p->plane1],
		      values, NULL, NULL, NULL);
	rtnVal = WlzHasIntersection(tmpObj1, tmpObj2, &errNum);
	WlzFreeObj(tmpObj1);
	WlzFreeObj(tmpObj2);
      }
      break;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnVal;
}
