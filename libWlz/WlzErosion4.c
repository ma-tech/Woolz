#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzErosion4_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzErosion4.c
* \author       Richard Baldock
* \date         September 2003
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
* \brief	Performs 4-connected erosion on domain objects.
* \ingroup	WlzMorphologyOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

static int intvsh_intv(WlzInterval *inta,
		       WlzInterval *intb,
		       WlzInterval *buff,
		       int *flag);
static int line_intsh_ints(WlzIntervalLine *inta,
			   WlzInterval *jntl,
			   int n,
			   WlzInterval *buff);
static int line_int_int(WlzIntervalLine *itvl,
			WlzIntervalLine *jtvl,
			WlzInterval *jp);
/* function:     WlzErosion4    */
/*! 
* \ingroup      WlzMorphologyOps
* \brief        4-connected erosion of a woolz domain object.2D objects
 only.This should not really be publicly accessible, but is present
 for historical reasons. Therefore the prototype does not appear in
 WlzProto.h. WlzErosion should be used to access 4-connected erosion.
*
* \return       Eroded object.
* \param    obj	Input object.
* \param    dstErr	Error return.
* \par      Source:
*                WlzErosion4.c
*/
WlzObject *WlzErosion4(
  WlzObject	*obj,
  WlzErrorNum	*dstErr)
{
  WlzObject		*erosobj=NULL;
  WlzDomain		domain;
  WlzValues		values;
  WlzIntervalDomain	*idmn1;
  WlzIntervalDomain	*idmn2;
  WlzIntervalLine	*itvl;
  WlzInterval		*jp, *ptr;
  WlzInterval		*jpw;
  WlzInterval		*buff = NULL;
  int			i, n, line;
  int			line1, lastln, totint;
  int			maxItvLn;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  values.core = NULL;
  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      /* check the domain */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      break;

    case WLZ_TRANS_OBJ:
    case WLZ_EMPTY_OBJ:
    case WLZ_3D_DOMAINOBJ:
    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    idmn2 = obj->domain.i;
    line1 = idmn2->line1 + 1;
    lastln = idmn2->lastln - 1;
    switch( idmn2->type ){

    case WLZ_INTERVALDOMAIN_INTVL:
      break;

    case WLZ_INTERVALDOMAIN_RECT:
      if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					   line1, lastln,
					   idmn2->kol1 + 1,
					   idmn2->lastkl - 1,
					   &errNum)) != NULL){
	return WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }

  /*
   * if the object is too thin after erosion the object will disappear
   */
  domain.core = NULL;
  values.core = NULL;
  if( errNum == WLZ_ERR_NONE ){
    if( (lastln < line1) || ((idmn2->lastkl - idmn2->kol1) < 2) ){
      return WlzMakeEmpty(dstErr);
    }
  }

  /*
   * reserve space for erosion object
   */
  if( errNum == WLZ_ERR_NONE ){
    totint =  2 * WlzIntervalCount(idmn2, NULL);
    if( (jp = (WlzInterval *) AlcMalloc(totint*sizeof(WlzInterval))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      jpw = jp;
    }
  }

  /* Make buffers with room for the maximum number of intervals in any line. */
  if( errNum == WLZ_ERR_NONE ){
    maxItvLn = WlzIDomMaxItvLn(idmn2);
    if(maxItvLn < 1){
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if((buff = (WlzInterval *)
                    AlcMalloc(sizeof(WlzInterval) * maxItvLn)) == NULL) {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if((errNum == WLZ_ERR_NONE) && 
     (idmn1 = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
				    line1, lastln,
				    idmn2->kol1+1, idmn2->lastkl-1,
				    &errNum)) ){
    idmn1->freeptr = AlcFreeStackPush(idmn1->freeptr, (void *)jp, NULL); 
    itvl = &(idmn2->intvlines[line1 - idmn2->line1]);
    ptr = buff;
    for(line = line1; line <= lastln; line++){
      n = line_int_int(itvl-1, itvl+1, ptr);
      i = line_intsh_ints(itvl, ptr, n, jp);
      WlzMakeInterval(line, idmn1, i, jp);
      jp += i;
      itvl++;
    }
    domain.i = idmn1;

    /* no used intervals implies empty object */
    if(jp == jpw){
      WlzFreeDomain(domain);
      return WlzMakeEmpty(dstErr);
    }
    erosobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			  NULL, NULL, &errNum);
  }
  else if( jp ){
    /* reaching here implies an error and jp != NULL */
    AlcFree( jp );
  }

  AlcFree(buff);

  if( dstErr ){
    *dstErr = errNum;
  }
  return erosobj;
}

/*
 * the same as intv_intv() in erosion exept that
 * this function shrink the interval of inta
 */
static int intvsh_intv(WlzInterval *inta,
		       WlzInterval *intb,
		       WlzInterval *buff,
		       int *flag)
{
  int r, left, right;

  r = inta->iright - 2;
  left = intb->ileft - 1;
  right = intb->iright - 1;
  if(r < right){
    *flag = 0;
    if(r < left)
      return(0);
    else {
      buff->ileft = WLZ_MAX(inta->ileft, left);
      if(buff->ileft > r)
	return(0);
      buff->iright = r;
    }
  } else {
    *flag = 1;
    if(right < inta->ileft)
      return(0);
    else {
      buff->ileft = WLZ_MAX(inta->ileft, left);
      if(buff->ileft > right)
	return(0);
      buff->iright = right;
    }
  }
  return(1);
}

/*
 * the same as line_int_int() in erosion exept that
 * the intervals of inta are shriked
 */
static int line_intsh_ints(WlzIntervalLine *inta,
			   WlzInterval *jntl,
			   int n,
			   WlzInterval *buff)
{
  WlzInterval *intl;
  WlzInterval *ptr = buff;
  int m;
  int t;

  /*
   * compare one interval by one interval
   */
  m = inta->nintvs;
  intl = inta->intvs;
  while(m != 0 && n != 0){
    if(intvsh_intv(intl, jntl, buff, &t))		
      buff++;
    if(t == 0){
      intl++;
      m--;
    } else {
      jntl++;
      n--;
    }
  }
  return((int)(buff-ptr));
}

static int line_int_int(WlzIntervalLine *itvl,
			WlzIntervalLine *jtvl,
			WlzInterval *jp)
{
  WlzInterval *intl, *jntl;
  WlzInterval *ptr = jp;
  int m, n;

  intl = itvl->intvs;
  m = itvl->nintvs;
  jntl = jtvl->intvs;
  n = jtvl->nintvs;
  while(m > 0 && n > 0){
    jp->ileft = WLZ_MAX(intl->ileft, jntl->ileft);
    jp->iright = WLZ_MIN(intl->iright, jntl->iright);
    if(jp->iright >= jp->ileft)
      jp++;
    if(intl->iright == jntl->iright){
      intl++;
      m--;
      jntl++;
      n--;
    } else if(intl->iright > jntl->iright){
      jntl++;
      n--;
    } else if(intl->iright < jntl->iright){
      intl++;
      m--;
    }
  }
  return((int)(jp - ptr));
}


