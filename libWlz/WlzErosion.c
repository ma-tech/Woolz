#pragma ident "MRC HGU $Id$"
/*!
* \file         libWlz/WlzErosion.c
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
* \brief	Morphological erosion of domain objects.
* \ingroup	WlzMorphologyOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

#define MAXLNITV 100

static int 			line_int_int(
				  WlzIntervalLine *inta,
				  WlzIntervalLine *intb,
				  WlzInterval *buff);
static int 			line_intv_arr(
				  WlzIntervalLine *itvl,
			 	  WlzInterval *buff,
			 	  int n,
			 	  WlzInterval *tmp);
static int 			intv_intv(WlzInterval *inta,
		     		  WlzInterval *intb,
		     		  WlzInterval *buff);
static int 			intv_arr(WlzInterval *intl,
		    		  WlzInterval *buff,
		    		  WlzInterval *tmp);

extern WlzObject 		*WlzErosion4(WlzObject *obj,
			      	  WlzErrorNum *wlzErr);

static WlzObject 		*WlzErosion3d(WlzObject *obj,
			          WlzConnectType  connectivity,
			       	  WlzErrorNum *wlzErr);

/*!
* \return	Eroded object or without values or NULL on error.
* \ingroup 	WlzMorphologyOps
* \brief	Calculates the morphological erosion of a woolz object
*		with a structuring element defined by the connectivity.
* \param	obj			Object to be eroded, must be a 2D
*					or 3D domain object (including a
*					WlzTransObj).
* \param	connectivity		Type of connectivity.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzErosion(
  WlzObject		*obj,
  WlzConnectType 	connectivity,
  WlzErrorNum		*dstErr)
{
  WlzObject 		*erosobj=NULL;
  WlzDomain 		domain;
  WlzValues		erosvalues,
  			values;
  WlzIntervalDomain 	*idmn;
  WlzIntervalLine 	*itvl;
  WlzInterval 		*jp, *jwp;
  WlzInterval 		buff[MAXLNITV], tmp[MAXLNITV];
  int			i, j;
  int 			inttot, *nitv;
  int 			line1, lastln, line;
  int 			k1, kol1,lastkl;
  int 			m, n, odd;
  int 			length, end;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  domain.core = NULL;
  values.core = NULL;

  /* check connectivity */
  switch( connectivity ){
  case WLZ_8_CONNECTED:
  case WLZ_4_CONNECTED:
  case WLZ_6_CONNECTED:
  case WLZ_18_CONNECTED:
  case WLZ_26_CONNECTED:
    break;

  default:
    errNum = WLZ_ERR_CONNECTIVITY_TYPE;
  }

  if( errNum == WLZ_ERR_NONE ){
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      /* check the domain */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      if( (connectivity == WLZ_4_CONNECTED) ||
	 (connectivity == WLZ_6_CONNECTED) ){
	return WlzErosion4(obj, dstErr);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzErosion3d(obj, connectivity, dstErr);

    case WLZ_TRANS_OBJ:
      if( (erosobj =
	   WlzErosion(obj->values.obj, connectivity, &errNum)) == NULL ){
	break;
      }
      erosvalues.obj = erosobj;
      return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, erosvalues,
			 NULL, NULL, dstErr);

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    idmn = obj->domain.i;
    switch( idmn->type ){

    case WLZ_INTERVALDOMAIN_INTVL:
      kol1 = idmn->kol1;
      lastkl = idmn->lastkl;
      line1 = idmn->line1;
      lastln = idmn->lastln;
      i = lastln - line1;
      if(i < 2 || lastkl - kol1 < 2){
	return WlzMakeMain(WLZ_EMPTY_OBJ, domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_INTERVALDOMAIN_RECT:
      kol1 = idmn->kol1 + 1;
      lastkl = idmn->lastkl - 1;
      line1 = idmn->line1 + 1;
      lastln = idmn->lastln - 1;
      if( kol1 > lastkl || line1 > lastln ){
	return WlzMakeEmpty(dstErr);
      }
      if( domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					    line1, lastln, kol1, lastkl,
					    &errNum) ){
	return WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_DOMAIN:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }

  /*
   * Otherwise, a interval-domain object - better test, then
   * reserve space for erosion object
   */
  if( errNum == WLZ_ERR_NONE ){
    if( (nitv = (int *) AlcMalloc(sizeof(int) * (i+1))) == NULL ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    nitv[i] = 0;
    inttot = WlzIntervalCount(idmn, &errNum);
    inttot += 2 * abs (inttot - i);
  }
  if( errNum == WLZ_ERR_NONE ){
    if( (jp = (WlzInterval *) AlcMalloc(inttot * sizeof(WlzInterval)))
       == NULL ){
      AlcFree((void *) nitv);
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      jwp = jp;
      odd = (i+1) & 01;
    }
  }

  /*
   * initialize
   * ----------
   */
  if( errNum == WLZ_ERR_NONE ){
    line1++;
    lastln--;
    k1 = kol1;
    kol1 = lastkl;
    lastkl = 0;
    itvl = idmn->intvlines;
    m = 0;
    itvl++;
    if(odd == 1)
      lastln--;
    for(line = line1; line <= lastln; line += 2){
      nitv[m] = 0;
      length = line_int_int(itvl, itvl+1, buff);
      for(j = -1; j <= 2; j += 3){
	end = line_intv_arr(itvl+j,buff,length,tmp);
	nitv[m] = 0;
	for(i = 0; i < end; i++){
	  if(tmp[i].iright - tmp[i].ileft > 1){
	    jp->ileft = tmp[i].ileft;
	    jp->iright = tmp[i].iright;
	    kol1 = WLZ_MIN(kol1, jp->ileft);
	    lastkl = WLZ_MAX(lastkl, jp->iright);
	    jp++;
	    nitv[m] += 1;
	  }
	}
	m++;
      }
      itvl += 2;
    }
    if(odd == 1){
      j = line_int_int(itvl-1, itvl, buff);
      length = line_intv_arr(itvl+1,buff,j,tmp);
      nitv[m] = 0;
      for(i = 0; i < length; i++){
	if(tmp[i].iright - tmp[i].ileft > 1){
	  jp->ileft = tmp[i].ileft;
	  jp->iright = tmp[i].iright;
	  kol1 = WLZ_MIN(kol1, jp->ileft);
	  lastkl = WLZ_MAX(lastkl, jp->iright);
	  jp++;
	  nitv[m] += 1;
	}
      }
      if(nitv[m] > 0)
	m++;
    }
  }
  /*
   * get first and last line
   */
  /*  - this tests the number of intervals  - zero implies empty */
  if( errNum == WLZ_ERR_NONE ){
    if(jp == jwp){
      AlcFree((char *) nitv);
      AlcFree((char *) jwp);
      return WlzMakeMain(WLZ_EMPTY_OBJ, domain, values, NULL, NULL, dstErr);
    }
    i = m - 1;
    while(nitv[i] == 0)
      i--;
    lastln = line1 + i;
    i = 0;
    while(nitv[i] == 0)
      i++;
    line1 += i;

    /* test for no lines */
    if(lastln < line1){
      AlcFree((char *) nitv);
      AlcFree((char *) jwp);
      return WlzMakeMain(WLZ_EMPTY_OBJ, domain, values, NULL, NULL, dstErr);
    }
  }

  /* make the new object */
  if( errNum == WLZ_ERR_NONE ){
    if( (domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					  line1, lastln, kol1+k1+1,
					  lastkl+k1-1, &errNum)) == NULL ){
      AlcFree((char *) nitv);
      AlcFree((char *) jwp);
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    if( (erosobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			       NULL, NULL, &errNum)) == NULL ){
      AlcFree((char *) nitv);
      AlcFree((char *) jwp);
      WlzFreeDomain(domain);
    }
  }
  
  if( errNum == WLZ_ERR_NONE ){
    jp = jwp;
    domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)jp, NULL);
    k1 = kol1;
    kol1 += 2;
    for(line = line1; line <= lastln; line++){
      jwp = jp;
      n = nitv[i];
      while (n-- > 0){
	jp->ileft -= k1;
	jp->iright -= kol1;
	jp++;
      }
      WlzMakeInterval(line, domain.i, nitv[i++], jwp);
    }
    AlcFree((void *) nitv);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return erosobj;
}

/*!
* \return	Number of intervals in the intersection.
* \ingroup	WlzMorphologyOps
* \brief	Finds the number of intersection intervals between two lines.
* \param	inta			Set of intervals on one line.
* \param	intb			Set of intervals on other line.
* \param	buff			Inverval buffer.
*/
static int line_int_int(WlzIntervalLine *inta,
			WlzIntervalLine *intb,
			WlzInterval *buff)
{
  WlzInterval *intl, *jntl;
  int i, m, n, t;

  /*
   * compare one interval by one interval
   */
  m = inta->nintvs;
  intl = inta->intvs;
  n = intb->nintvs;
  jntl = intb->intvs;
  i = 0;
  while(m != 0 && n != 0){
    t = intv_intv(intl, jntl, buff);		
    if(t == -1){
      intl++;
      m--;
    } else if (t == 2){
      jntl++;
      n--;
    } else {
      if(t == 0){
	intl++;
	m--;
      } else {
	jntl++;
	n--;
      }
      buff++;
      i++;
    }
  }
  return(i);
}

/*!
* \return	Number of intersection intervals.
* \ingroup	WlzMorphologyOps
* \brief	Finds the intersection intervals between two lines.
*		One set of intervals is pointed out by itvl the other n
*		intervals is in array left and right length is the result
*		number of intersection intervals in array l and r.
* \param	itvl
* \param	buff
* \param	n
* \param	tmp
*/
static int line_intv_arr(WlzIntervalLine *itvl,
			 WlzInterval *buff,
			 int n,
			 WlzInterval *tmp)
{
  WlzInterval *intl;
  int i, j, m, t;

  /*
   * compare one interval by one interval
   */
  intl = itvl->intvs;
  m = itvl->nintvs;
  i = 0;
  j = 0;
  while(m != 0 && j != n){
    t = intv_arr(intl, &buff[j], &tmp[i]);		
    if(t == -1){
      intl++;
      m--;
    } else if (t == 2)
      j++;
    else {
      if(t == 0){
	intl++;
	m--;
      } else j++;
      i++;
    }
  }
  return(i);
}

/*!
* \return	<ul>
*		<li> -1: inta is in far left, there is no intersection point
* 	  	<li> 0: there is an intersection and intb at the right of inta
* 	  	<li> 1: there is an intersection and inta at the right of intb
* 		<li> 2: intb is in far right, there is no intersection point
*		</ul>
* \brief	Finds the intersection of two intervals.
* \param	inta			First interval.
* \param	intb			Second interval.
* \param	buff			Interval buffer.
*/
static int intv_intv(WlzInterval *inta,
		     WlzInterval *intb,
		     WlzInterval *buff)
{
  if(inta->iright < intb->ileft)
    return(-1);
  if(inta->ileft > intb->iright)
    return(2);
  buff->ileft = WLZ_MAX(inta->ileft, intb->ileft);
  if(inta->iright < intb->iright){
    buff->iright = inta->iright;
    return(0);
  } else {
    buff->iright = intb->iright;
    return(1);
  }
}
	
/*
 * find intersection of an interval and a section(left, right)
 * return -1: inta is in far left, there is no intersection point
 * return 2: left is in far right, there is no intersection point
 * 	  0: there is an intersection and intl is at the left of "right"
 * 	  1: there is an intersection and "right" is at the left of intl
 */
static int intv_arr(WlzInterval *intl,
		    WlzInterval *buff,
		    WlzInterval *tmp)
{
  if(intl->iright < buff->ileft)
    return(-1);
  if(intl->ileft > buff->iright)
    return(2);
  tmp->ileft = WLZ_MAX(intl->ileft, buff->ileft);
  if(intl->iright < buff->iright){
    tmp->iright = intl->iright;
    return(0);
  } else {
    tmp->iright = buff->iright;
    return(1);
  }
}

/*!
* \return	Eroded object or NULL on error.
* \ingroup	WlzMorphologyOps
* \brief	Erosion of 3D objects. Only called from WlzErosion().
* \param	obj			Object to be eroded.
* \param	connectivity		Required connectivity.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzErosion3d(
  WlzObject 		*obj,
  WlzConnectType 	connectivity,
  WlzErrorNum		*dstErr)
{
  WlzObject	*new_obj=NULL;
  WlzObject	*start_obj[3], *dest_obj[3], *tmp_obj;
  WlzDomain	domain;
  WlzValues	values;
  int		p, nplanes;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* no need to check object pointer or type
     but do need to check the 3D domain type */
  switch( obj->domain.p->type ){

  case WLZ_PLANEDOMAIN_DOMAIN:
    nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 1;
    break;

  default:
  case WLZ_PLANEDOMAIN_BOUNDLIST:
    errNum = WLZ_ERR_DOMAIN_TYPE;
    break;
  }

  /* create a new 3D object to hold the eroded domains
     use WlzStandardPlaneDomain to trim leading and trailing
     empty domains */
  if( errNum == WLZ_ERR_NONE ){
    if( domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				      obj->domain.p->plane1,
				      obj->domain.p->lastpl,
				      obj->domain.p->line1,
				      obj->domain.p->lastln,
				      obj->domain.p->kol1,
				      obj->domain.p->lastkl, &errNum) ){
      for(p=0; p < 3; p++){
	domain.p->voxel_size[p] = obj->domain.p->voxel_size[p];
      }
    }
  }
  if( errNum == WLZ_ERR_NONE ){
    values.core = NULL;
    new_obj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			  NULL, NULL, &errNum);
  }

  /* foreach plane erode as required by connectivity */
  if( errNum == WLZ_ERR_NONE ){
    domain.core = NULL;
    values.core = NULL;
    start_obj[0] = WlzAssignObject(WlzMakeEmpty(NULL), NULL);
    start_obj[1] = WlzAssignObject(WlzMakeEmpty(NULL), NULL);
    start_obj[2] =
      WlzAssignObject (WlzMakeMain(WLZ_2D_DOMAINOBJ,
				   obj->domain.p->domains[0],
				   values, NULL, NULL, NULL), NULL);
    for(p=0; p < nplanes; p++){
      WlzFreeObj( start_obj[0] );
      start_obj[0] = start_obj[1];
      start_obj[1] = start_obj[2];
      if( (p < (nplanes-1)) && obj->domain.p->domains[p+1].core ){
	start_obj[2] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				   obj->domain.p->domains[p+1],
				   values, NULL, NULL, NULL);
      }
      else {
	start_obj[2] = WlzMakeEmpty(NULL);
      }
      start_obj[2] = WlzAssignObject(start_obj[2], NULL);

      switch( connectivity ){

      case WLZ_8_CONNECTED:
      case WLZ_4_CONNECTED:
	if( dest_obj[1] = WlzErosion(start_obj[1], connectivity, NULL) ){
	  if( dest_obj[1]->type == WLZ_EMPTY_OBJ ){
	    new_obj->domain.p->domains[p].core = NULL;
	  }
	  else {
	    new_obj->domain.p->domains[p] =	
	      WlzAssignDomain(dest_obj[1]->domain, NULL);
	  }
	  WlzFreeObj( dest_obj[1] );
	}
	else {
	  new_obj->domain.p->domains[p].core = NULL;
	}
	break;

      case WLZ_6_CONNECTED:
	dest_obj[0] = WlzAssignObject(start_obj[0], NULL);
	dest_obj[1] = WlzErosion(start_obj[1], WLZ_4_CONNECTED, NULL);
	dest_obj[2] = WlzAssignObject(start_obj[2], NULL);
	if( tmp_obj = WlzIntersectN(3, dest_obj, 0, NULL) ){
	  if( tmp_obj->type == WLZ_EMPTY_OBJ){
	    new_obj->domain.p->domains[p].core = NULL;
	  }
	  else{
	    new_obj->domain.p->domains[p] = WlzAssignDomain(tmp_obj->domain,
							    NULL);
	  }
	  WlzFreeObj( tmp_obj );
	}
	else {
	  new_obj->domain.p->domains[p].core = NULL;
	}
	WlzFreeObj( dest_obj[0] );
	WlzFreeObj( dest_obj[1] );
	WlzFreeObj( dest_obj[2] );
	break;

      case WLZ_18_CONNECTED:
	dest_obj[0] = WlzErosion(start_obj[0], WLZ_4_CONNECTED, NULL);
	dest_obj[1] = WlzErosion(start_obj[1], WLZ_8_CONNECTED, NULL);
	dest_obj[2] = WlzErosion(start_obj[2], WLZ_4_CONNECTED, NULL);
	if( tmp_obj = WlzIntersectN(3, dest_obj, 0, NULL) ){
	  if( tmp_obj->type == WLZ_EMPTY_OBJ){
	    new_obj->domain.p->domains[p].core = NULL;
	  }
	  else{
	    new_obj->domain.p->domains[p] = WlzAssignDomain(tmp_obj->domain,
							    NULL);
	  }
	  WlzFreeObj( tmp_obj );
	}
	else {
	  new_obj->domain.p->domains[p].core = NULL;
	}
	WlzFreeObj( dest_obj[0] );
	WlzFreeObj( dest_obj[1] );
	WlzFreeObj( dest_obj[2] );
	break;

      case WLZ_26_CONNECTED:
	dest_obj[0] = WlzErosion(start_obj[0], WLZ_8_CONNECTED, NULL);
	dest_obj[1] = WlzErosion(start_obj[1], WLZ_8_CONNECTED, NULL);
	dest_obj[2] = WlzErosion(start_obj[2], WLZ_8_CONNECTED, NULL);
	if( tmp_obj = WlzIntersectN(3, dest_obj, 0, NULL) ){
	  if( tmp_obj->type == WLZ_EMPTY_OBJ){
	    new_obj->domain.p->domains[p].core = NULL;
	  }
	  else{
	    new_obj->domain.p->domains[p] = WlzAssignDomain(tmp_obj->domain,
							    NULL);
	  }
	  WlzFreeObj( tmp_obj );
	}
	else {
	  new_obj->domain.p->domains[p].core = NULL;
	}
	WlzFreeObj( dest_obj[0] );
	WlzFreeObj( dest_obj[1] );
	WlzFreeObj( dest_obj[2] );
	break;

      }
    }
    WlzFreeObj( start_obj[0] );
    WlzFreeObj( start_obj[1] );
    WlzFreeObj( start_obj[2] );

    WlzStandardPlaneDomain(new_obj->domain.p, NULL);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return new_obj;
}
