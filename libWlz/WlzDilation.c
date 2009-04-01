#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDilation_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzDilation.c
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
* \brief	Functions for dilating objects with spatial domains.
* \ingroup	WlzMorphologyOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>

static int line_int_dil(WlzIntervalLine *inta,
			WlzIntervalLine *intb,
			WlzInterval 	*cc);

static int line_int_cp(WlzIntervalLine *inta,
		       WlzIntervalLine *intb,
		       WlzInterval 	*cc);

static int set_line_dil(WlzIntervalLine *itvl,
			WlzInterval 	*jntl);

static int set_line_cp(WlzIntervalLine *itvl,
			WlzInterval 	*jntl);

static int line_arr_dil(WlzIntervalLine *inta,
			WlzInterval 	*jntl,
			int 		n,
			WlzInterval 	*cc);

static WlzObject *WlzDilation3d(WlzObject *obj,
				WlzConnectType 	connectivity,
				WlzErrorNum	*dstErr);
static WlzObject *WlzDilation4(WlzObject *obj,
			       WlzErrorNum	*dstErr);
/*!
* \return	Dilated object.
* \ingroup	WlzMorphologyOps
* \brief	Dilate the given object using the given connectivity type.
* 		Since the dilated object is bigger than the original, the
*		size of the valuetable may be smaller than the dilated object.
*		User has to take fully responsibility for using grey value of
*		dilated object.
* \param	obj			Given object.
* \param	connectivity		Required type of conectivity.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzDilation(
  WlzObject 		*obj,
  WlzConnectType 	connectivity,
  WlzErrorNum		*dstErr)
{
  WlzObject 		*dilatobj=NULL;
  WlzDomain		domain;
  WlzValues		dilatvalues;
  WlzIntervalDomain 	*idmn;
  WlzIntervalLine 	*itvl;
  WlzInterval 		*buff = NULL;
  WlzInterval 		*jp;
  int 			i, j;
  int			maxItvLn;
  int 			inttot, nitv;
  int 			line1, lastln, line;
  int 			kol1,lastkl;
  int 			odd, length;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* Check object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    domain.core = NULL;
    dilatvalues.core = NULL;
    switch( obj->type ){

    case WLZ_2D_DOMAINOBJ:
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
	break;
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzDilation3d(obj, connectivity, dstErr);

    case WLZ_TRANS_OBJ:
      if( (dilatobj =
	   WlzDilation(obj->values.obj, connectivity, &errNum)) == NULL ){
	break;
      }
      dilatvalues.obj = dilatobj;
      return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, dilatvalues,
			 NULL, NULL, dstErr);

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* check connectivity */
  if( errNum == WLZ_ERR_NONE ){
    switch( connectivity ){
    case WLZ_8_CONNECTED:
    case WLZ_18_CONNECTED:
    case WLZ_26_CONNECTED:
      connectivity = WLZ_8_CONNECTED;
      break;

    case WLZ_4_CONNECTED:
    case WLZ_6_CONNECTED:
      connectivity = WLZ_4_CONNECTED;
      return WlzDilation4(obj, dstErr);
      
    default:
      errNum = WLZ_ERR_CONNECTIVITY_TYPE;
    }
  }

  /* check the domain */
  if( errNum == WLZ_ERR_NONE ){
    idmn = obj->domain.i;
    switch( idmn->type ){

    case WLZ_INTERVALDOMAIN_INTVL:
      kol1 = idmn->kol1;
      lastkl = idmn->lastkl;
      line1 = idmn->line1;
      lastln = idmn->lastln;	
      break;

    case WLZ_INTERVALDOMAIN_RECT:
      kol1 = idmn->kol1 - 1;
      lastkl = idmn->lastkl + 1;
      line1 = idmn->line1 - 1;
      lastln = idmn->lastln + 1;
      if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_RECT,
					    line1, lastln, kol1, lastkl,
					    &errNum)) != NULL){
	return WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, dilatvalues,
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

  if( errNum == WLZ_ERR_NONE ){
    i = lastln - line1 + 3;
    inttot = WlzIntervalCount(idmn, &errNum);
    /* I don't understand this interval count estimate - it must fail
       occasionally but the #new intervals <= 3 * #old intervals therefore
       use this formula and shrink the storage requirement after the
       calculation using realloc if it seems important */
    /* Was inttot += 10 * MAXLNITV, with MAXLNITV = 300 */
    inttot *= 3;
    if( errNum == WLZ_ERR_NONE ){
      if((jp = (WlzInterval *)AlcMalloc(inttot *
                                        sizeof(WlzInterval))) != NULL)
      {
	odd = i & 01;
	if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					     line1-1, lastln+1,
					     kol1-1, lastkl+1,
					     &errNum)) != NULL){
	  domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)jp,
	  				       NULL);
	  if((dilatobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, dilatvalues,
				     NULL, NULL, &errNum)) == NULL){
	    WlzFreeDomain(domain);
	  }
	}
      }
      else {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }

  /* Make buffer with room for the maximum number of intervals in any line
   * x 3. This replaces a previous hard coded maximum number (300), but
   * I'm not convinced that this is right. */
  if( errNum == WLZ_ERR_NONE ){
    maxItvLn = WlzIDomMaxItvLn(domain.i) * 3;
    if(maxItvLn < 1){
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if((buff = (WlzInterval *)
                    AlcMalloc(sizeof(WlzInterval) *maxItvLn)) == NULL) {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    /*
     * initialize
     * ----------
     */
    itvl = idmn->intvlines;
    if(odd == 1){
      lastln--;
    }
    /*
     * get first two lines because they only depend on one or two lines
     */
    line = line1 - 1;
    nitv = set_line_dil(itvl, jp);
    WlzMakeInterval(line++, domain.i, nitv, jp);
    jp += nitv;
    /*
     * if there is only one line make life easier
     */
    if(idmn->line1 == idmn->lastln){
      for(i = 0; i < 2; i++){
	nitv = set_line_dil(itvl, jp);
	WlzMakeInterval(line++, domain.i, nitv, jp);
	jp += nitv;
      }
    }
    else{
      nitv = line_int_dil(itvl, itvl+1, jp);
      WlzMakeInterval(line++, domain.i, nitv, jp);
      jp += nitv;
      itvl++;
      /*
       * get the rest lines which need combine three lines
       */
      for(line = line1 + 1; line < lastln; ){
	length = line_int_dil(itvl, itvl+1, buff);
	for(j = -1; j <= 2; j += 3){
	  nitv = line_arr_dil(itvl+j, buff, length, jp);
	  WlzMakeInterval(line++, domain.i, nitv, jp);
	  jp += nitv;
	}
	itvl += 2;
      }
      line = lastln;
      if(odd == 1){
	j = line_int_dil(itvl-1, itvl, buff);
	nitv = line_arr_dil(itvl+1, buff, j, jp);
	WlzMakeInterval(line++, domain.i, nitv, jp);
	jp += nitv;
      }
      /*
       * get final two lines
       */
      if(odd == 0)
	itvl--;
      nitv = line_int_dil(itvl, itvl+1, jp);
      WlzMakeInterval(line++, domain.i, nitv, jp);
      jp += nitv;
      itvl++;
      nitv = set_line_dil(itvl, jp);
      WlzMakeInterval(line, domain.i, nitv, jp);
      jp += nitv;
    }
  }

  AlcFree(buff);

  if( dstErr ){
    *dstErr = errNum;
  }
  return(dilatobj);
}

/*!
* \return	Number of intersection intervals in interval array cc.
* \ingroup	WlzMorphologyOps
* \brief	Finds the union intervals between two lines. One set of
* 		intervals is pointed out by intva the other intervals is
* 		pointed out by intvb length is the result number of
* 		intersection intervals in interval array cc.
* \param	inta			First set of intervals.
* \param	intb			Second set of intervals.
* \param	cc			Intersection intervals.
*/
static int line_int_dil(WlzIntervalLine 	*inta,
			WlzIntervalLine 	*intb,
			WlzInterval 	*cc)
{
  WlzInterval *intl, *jntl;
  WlzInterval *tmp = cc;
  int m, n;
  int i;

  intl = inta->intvs;
  m = inta->nintvs;
  jntl = intb->intvs;
  n = intb->nintvs;
  if(m == 0)
    return(set_line_dil(intb, cc));
  else if(n == 0)
    return(set_line_dil(inta, cc));
  if(intl->ileft > jntl->iright){
    cc->ileft = jntl->ileft;
    cc->iright = jntl->iright + 2;
    jntl++;
    n--;
  } else if(intl->iright < jntl->ileft){
    cc->ileft = intl->ileft;
    cc->iright = intl->iright + 2;
    intl++;
    m--;
  } else {
    cc->ileft = WLZ_MIN(intl->ileft, jntl->ileft);
    cc->iright = WLZ_MAX(intl->iright, jntl->iright) + 2;
    intl++;
    jntl++;
    m--;
    n--;
  }
  while(m > 0 && n > 0){
    if(intl->ileft < jntl->ileft){
      if(cc->iright + 1 >= intl->ileft)
	cc->iright = WLZ_MAX(cc->iright,intl->iright+2);
      else {
	cc++;
	cc->iright = intl->iright + 2;
	cc->ileft = intl->ileft;
      }
      intl++;
      m--;
    } else if(intl->ileft > jntl->ileft){
      if(cc->iright + 1 >= jntl->ileft)
	cc->iright = WLZ_MAX(cc->iright,jntl->iright+2);
      else {
	cc++;
	cc->iright = jntl->iright + 2;
	cc->ileft = jntl->ileft;
      }
      jntl++;
      n--;
    } else {
      i = WLZ_MAX(intl->iright, jntl->iright) + 2;
      if(cc->iright + 1 >= intl->ileft){
	cc->iright = i;
      } else {
	cc++;
	cc->ileft = intl->ileft;
	cc->iright = i;
      }
      intl++;
      jntl++;
      m--;
      n--;
    }
  }
  if(n > 0){
    m = n;
    intl = jntl;
  }
  while(m-- > 0){ if(cc->iright + 1 < intl->ileft){
    cc++;
    cc->ileft = intl->ileft;
    cc->iright = intl->iright + 2;
  } else {
    cc->iright = WLZ_MAX(cc->iright, intl->iright + 2);
  }
  intl++;
  }
  cc++;
  return((int)(cc-tmp));
}

/*!
* \return	The number of intersection intervals in interval array cc.
* \ingroup	WlzMorphologyOps
* \brief	Finds the union intervals between two lines. One set of
*		intervals is pointed out by intva the other intervals is
*		pointed out by intvb. The length is the result number of
*		intersection intervals in interval array cc - non dilate
*		version for 4-connected dilation.
* \param	inta			First set of intervals.
* \param	intb			Second set of intervals.
* \param	cc			Intersection intervals.
*/
static int line_int_cp(WlzIntervalLine 	*inta,
		       WlzIntervalLine 	*intb,
		       WlzInterval 	*cc)
{
  WlzInterval *intl, *jntl;
  WlzInterval *tmp = cc;
  int m, n;
  int i;

  intl = inta->intvs;
  m = inta->nintvs;
  jntl = intb->intvs;
  n = intb->nintvs;
  if(m == 0)
    return(set_line_cp(intb, cc));
  else if(n == 0)
    return(set_line_cp(inta, cc));
  if(intl->ileft > jntl->iright){
    cc->ileft = jntl->ileft + 1;
    cc->iright = jntl->iright + 1;
    jntl++;
    n--;
  } else if(intl->iright < jntl->ileft){
    cc->ileft = intl->ileft + 1;
    cc->iright = intl->iright + 1;
    intl++;
    m--;
  } else {
    cc->ileft = WLZ_MIN(intl->ileft, jntl->ileft) + 1;
    cc->iright = WLZ_MAX(intl->iright, jntl->iright) + 1;
    intl++;
    jntl++;
    m--;
    n--;
  }
  while(m > 0 && n > 0){
    if(intl->ileft < jntl->ileft){
      if(cc->iright >= intl->ileft)
	cc->iright = WLZ_MAX(cc->iright,intl->iright+1);
      else {
	cc++;
	cc->iright = intl->iright + 1;
	cc->ileft = intl->ileft + 1;
      }
      intl++;
      m--;
    } else if(intl->ileft > jntl->ileft){
      if(cc->iright >= jntl->ileft)
	cc->iright = WLZ_MAX(cc->iright,jntl->iright+1);
      else {
	cc++;
	cc->iright = jntl->iright + 1;
	cc->ileft = jntl->ileft + 1;
      }
      jntl++;
      n--;
    } else {
      i = WLZ_MAX(intl->iright, jntl->iright) + 1;
      if(cc->iright >= intl->ileft){
	cc->iright = i;
      } else {
	cc++;
	cc->ileft = intl->ileft + 1;
	cc->iright = i;
      }
      intl++;
      jntl++;
      m--;
      n--;
    }
  }
  if(n > 0){
    m = n;
    intl = jntl;
  }
  while(m-- > 0){ if(cc->iright < intl->ileft){
    cc++;
    cc->ileft = intl->ileft + 1;
    cc->iright = intl->iright + 1;
  } else {
    cc->iright = WLZ_MAX(cc->iright, intl->iright + 1);
  }
  intl++;
  }
  cc++;
  return((int)(cc-tmp));
}

/*!
* \return	Number of intervals.
* \ingroup	WlzMorphologyOps
* \brief	Put dilated intervals of itvl into interval array jntl.
* \param	itvl			Source of dilated intervals.
* \param	jntl			Destination for dilated intervals.
*/
static int set_line_dil(WlzIntervalLine 	*itvl,
			WlzInterval 	*jntl)
{
  WlzInterval *intl, *tmp = jntl;
  int m;

  /*
   * compare one interval by one interval
   */
  intl = itvl->intvs;
  m = itvl->nintvs;
  if(m == 0)
    return(0);
  jntl->ileft = intl->ileft;
  jntl->iright = intl->iright + 2;
  if(m == 1)
    return(1);
  intl++;
  m--;
  while(m-- > 0){
    if(jntl->iright + 1 < intl->ileft){
      jntl++;
      jntl->ileft = intl->ileft;
      jntl->iright = intl->iright + 2;
    } else {
      jntl->iright = intl->iright + 2;
    }
    intl++;
  }
  jntl++;
  return((int) (jntl-tmp));
}

/*!
* \return	Number of intervals.
* \ingroup	WlzMorphologyOps
* \brief	Copy intervals of itvl into interval array jntl - for
* 		4-connected dilation.
* \param	itvl			Source of intervals.
* \param	jntl			Destination of intervals.
*/
static int set_line_cp(WlzIntervalLine 	*itvl,
		       WlzInterval 	*jntl)
{
  WlzInterval	*intl;
  int		n;

  /*
   * copy interval array - incrementing values to account for the
   * new domain offset
   */
  n = itvl->nintvs;
  intl = itvl->intvs;
  while( n-- > 0 ){
    jntl->ileft = intl->ileft + 1;
    jntl->iright = intl->iright + 1;
    intl++;
    jntl++;
  }

  return itvl->nintvs;
}

/*!
* \return
* \ingroup	WlzMorphologyOps
* \brief	Put the union of intervals of jntl and the dilated intrvals of
* 		inta into cc. Note: inta has not been dilated but jntl has.
* \param	inta
* \param	jntl
* \param	n
* \param	cc
*/
static int line_arr_dil(WlzIntervalLine 	*inta,
			WlzInterval 	*jntl,
			int 			n,
			WlzInterval 	*cc)
{
  WlzInterval *intl;
  WlzInterval *tmp = cc;
  int i, m;

  intl = inta->intvs;
  m = inta->nintvs;
  if(n == 0)
    return(set_line_dil(inta, cc));
  else if(m == 0){
    i = n;
    while(n-- > 0){
      cc->ileft = jntl->ileft;
      cc->iright = jntl->iright;
      cc++;
      jntl++;
    }
    return(i);
  }
  if(intl->ileft > jntl->iright){
    cc->ileft = jntl->ileft;
    cc->iright = jntl->iright;
    jntl++;
    n--;
  } else if(intl->iright < jntl->ileft){
    cc->ileft = intl->ileft;
    cc->iright = intl->iright + 2;
    intl++;
    m--;
  } else {
    cc->ileft=(intl->ileft>jntl->ileft) ? jntl->ileft : intl->ileft;
    cc->iright = WLZ_MAX(intl->iright + 2, jntl->iright);
    intl++;
    jntl++;
    m--;
    n--;
  }
  while(m > 0 && n > 0){
    if(intl->ileft < jntl->ileft){
      if(cc->iright + 1 >= intl->ileft)
	cc->iright = WLZ_MAX(cc->iright,intl->iright+2);
      else {
	cc++;
	cc->iright = intl->iright + 2;
	cc->ileft = intl->ileft;
      }
      intl++;
      m--;
    } else if(intl->ileft > jntl->ileft){
      if(cc->iright + 1 >= jntl->ileft)
	cc->iright = WLZ_MAX(cc->iright,jntl->iright);
      else {
	cc++;
	cc->iright = jntl->iright;
	cc->ileft = jntl->ileft;
      }
      jntl++;
      n--;
    } else {
      i = WLZ_MAX(intl->iright+2, jntl->iright);
      if(cc->iright + 1 >= intl->ileft){
	cc->iright = i;
      } else {
	cc++;
	cc->ileft = intl->ileft;
	cc->iright = i;
      }
      intl++;
      jntl++;
      m--;
      n--;
    }
  }
  while(m-- > 0){ if(cc->iright + 1 < intl->ileft){
    cc++;
    cc->ileft = intl->ileft;
    cc->iright = intl->iright + 2;
  } else {
    cc->iright = WLZ_MAX(cc->iright, intl->iright + 2);
  }
  intl++;
  }
  while(n-- > 0){ if(cc->iright + 1 < jntl->ileft){
    cc++;
    cc->ileft = jntl->ileft;
    cc->iright = jntl->iright;
  } else {
    cc->iright = WLZ_MAX(cc->iright, jntl->iright);
  }
  jntl++;
  }
  cc++;
  return((int)(cc-tmp));
}

/*!
* \return
* \ingroup	WlzMorphologyOps
* \brief	Four connected dilation.
* \param	obj
* \param	dstErr
*/
static WlzObject *WlzDilation4(
  WlzObject *obj,
  WlzErrorNum	*dstErr)
{
  WlzObject	*dilatobj=NULL;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzDomain	domain;
  WlzValues	values;
  WlzInterval	*jp;
  int		maxItvLn;
  int		kol1, lastkl, line1, lastln;
  int		line, width, inttot, nitv;
  WlzIntervalDomain 	*idmn;
  WlzIntervalLine 	*itvl;
  WlzInterval 		*buff = NULL;
  int 			length;

  /* only need to check the domain type to set the dilated object */
  values.core = NULL;
  switch( obj->domain.i->type ){
  case WLZ_INTERVALDOMAIN_INTVL:
    idmn = obj->domain.i;
    kol1 = idmn->kol1;
    lastkl = idmn->lastkl;
    line1 = idmn->line1;
    lastln = idmn->lastln;	
    break;

  case WLZ_INTERVALDOMAIN_RECT:
    /* calculate explicitly for speed */
    kol1 = obj->domain.i->kol1;
    lastkl = obj->domain.i->lastkl;
    line1 = obj->domain.i->line1;
    lastln = obj->domain.i->lastln;
    width = lastkl - kol1 + 3;
    if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					     line1-1, lastln+1,
					     kol1-1, lastkl+1,
					     &errNum)) != NULL){
      inttot = lastln - line1 + 3;
      if((jp = (WlzInterval *) AlcMalloc(inttot *
                                         sizeof(WlzInterval))) != NULL){
	domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)jp,
					     NULL);
	/* first line */
	line = line1 - 1;
	jp->ileft = 1;
	jp->iright = width - 2;
	WlzMakeInterval(line++, domain.i, 1, jp++);
	/* middle lines */
	while( line <= lastln ){
	  jp->ileft = 0;
	  jp->iright = width - 1;
	  WlzMakeInterval(line++, domain.i, 1, jp++);
	}
	/* last line */
	jp->ileft = 1;
	jp->iright = width - 2;
	WlzMakeInterval(line++, domain.i, 1, jp++);
	return WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			   NULL, NULL, dstErr);
      }
      else {
	errNum = WLZ_ERR_MEM_ALLOC;
	WlzFreeDomain(domain);
      }
    }
    break;

  default:
    errNum = WLZ_ERR_DOMAIN_TYPE;
    break;
  }

  /* now do the interval domain case */
  if( errNum == WLZ_ERR_NONE ){
    inttot = WlzIntervalCount(idmn, &errNum);
    /* I don't understand this interval count estimate - it must fail
       occasionally but the #new intervals <= 3 * #old intervals therefore
       use this formula and shrink the storage requirement after the
       calculation using realloc if it seems important */
    /* Was inttot += 10 * MAXLNITV, with MAXLNITV = 300 */
    inttot *= 3;
    if( errNum == WLZ_ERR_NONE ){
      if((jp = (WlzInterval *) AlcMalloc(inttot *
                                         sizeof(WlzInterval))) != NULL)
      {
	if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					     line1-1, lastln+1,
					     kol1-1, lastkl+1,
					     &errNum)) != NULL){
	  domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)jp,
	  				       NULL);
	  if((dilatobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
				     NULL, NULL, &errNum)) == NULL){
	    WlzFreeDomain(domain);
	  }
	}
      }
      else {
	errNum = WLZ_ERR_MEM_ALLOC;
      }
    }
  }

  /* Make buffer with room for the maximum number of intervals in any line
   * x 3. This replaces a previous hard coded maximum number (300), but
   * I'm not convinced that this is right. */
  if( errNum == WLZ_ERR_NONE ){
    maxItvLn = WlzIDomMaxItvLn(domain.i) * 3;
    if(maxItvLn < 1){
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if((buff = (WlzInterval *)
                    AlcMalloc(sizeof(WlzInterval) * maxItvLn)) == NULL) {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    /*
     * initialize
     */
    itvl = idmn->intvlines;

    /*
     * get first two lines because they only depend on one or two lines
     */
    line = line1 - 1;
    nitv = set_line_cp(itvl, jp);
    WlzMakeInterval(line++, domain.i, nitv, jp);
    jp += nitv;
    /*
     * if there is only one line make life easier
     */
    if( idmn->line1 == idmn->lastln ){
      nitv = set_line_dil(itvl, jp);
      WlzMakeInterval(line++, domain.i, nitv, jp);
      jp += nitv;
      nitv = set_line_cp(itvl, jp);
      WlzMakeInterval(line++, domain.i, nitv, jp);
      jp += nitv;
    }
    else{
      /* first line of original domain */
      length = set_line_cp(itvl+1, buff);
      nitv = line_arr_dil(itvl, buff, length, jp);
      WlzMakeInterval(line++, domain.i, nitv, jp);
      jp += nitv;
      itvl++;

      /* all lines that have at least one following line */
      while( (lastln - line) ){
	length = line_int_cp(itvl-1, itvl+1, buff);
	nitv = line_arr_dil(itvl, buff, length, jp);
	WlzMakeInterval(line++, domain.i, nitv, jp);
	jp += nitv;
	itvl++;
      }

      /* final line of original domain */
      length = set_line_cp(itvl-1, buff);
      nitv = line_arr_dil(itvl, buff, length, jp);
      WlzMakeInterval(line++, domain.i, nitv, jp);
      jp += nitv;

      /* final line of the new domain */
      nitv = set_line_cp(itvl, jp);
      WlzMakeInterval(line, domain.i, nitv, jp);
    }
  }

  AlcFree(buff);

  if( dstErr ){
    *dstErr = errNum;
  }
  return dilatobj;
}

/*!
* \return	Dilated object.
* \ingroup	WlzMorphologyOps
* \brief	Dilation of a 3D object.
* \param	obj			Given object.
* \param	connectivity		Required connectivity.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzDilation3d(
  WlzObject 		*obj,
  WlzConnectType 	connectivity,
  WlzErrorNum		*dstErr)
{
  WlzObject		*dilatobj=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  WlzObject	*start_obj[3], *dest_obj[3], *tmp_obj;
  WlzDomain	domain;
  WlzValues	values;
  int		p, nplanes;

  /* no need to check object pointer or type
     but do need to check the 3D domain type */
  switch( obj->domain.p->type ){

  case WLZ_PLANEDOMAIN_DOMAIN:
    nplanes = obj->domain.p->lastpl - obj->domain.p->plane1 + 3;
    break;

  default:
  case WLZ_PLANEDOMAIN_BOUNDLIST:
    errNum = WLZ_ERR_DOMAIN_TYPE;
    break;
  }

  /* create a new 3D object to hold the dilated domains
     use WlzStandardPlaneDomain to trim leading and trailing
     empty domains */
  if(errNum == WLZ_ERR_NONE){
    if((domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				      obj->domain.p->plane1 - 1,
				      obj->domain.p->lastpl + 1,
				      obj->domain.p->line1 - 1,
				      obj->domain.p->lastln + 1,
				      obj->domain.p->kol1 - 1,
				      obj->domain.p->lastkl + 1,
				      &errNum)) != NULL){
      for(p=0; p < 3; p++){
	domain.p->voxel_size[p] = obj->domain.p->voxel_size[p];
      }
    }
  }
  if(errNum == WLZ_ERR_NONE){
    values.core = NULL;
    dilatobj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			  NULL, NULL, &errNum);
  }

  /* foreach plane dilate as required by connectivity */
  if(errNum == WLZ_ERR_NONE){
    domain.core = NULL;
    values.core = NULL;
    start_obj[0] = WlzAssignObject(WlzMakeEmpty(NULL), NULL);
    start_obj[1] = WlzAssignObject(WlzMakeEmpty(NULL), NULL);
    start_obj[2] = WlzAssignObject(WlzMakeEmpty(NULL), NULL);

    for(p=0; p < nplanes; p++){
      WlzFreeObj( start_obj[0] );
      start_obj[0] = start_obj[1];
      start_obj[1] = start_obj[2];
      if( (p < (nplanes-2)) && obj->domain.p->domains[p].core ){
	start_obj[2] = WlzMakeMain(WLZ_2D_DOMAINOBJ,
				   obj->domain.p->domains[p],
				   values, NULL, NULL, NULL);
      }
      else {
	start_obj[2] = WlzMakeEmpty(NULL);
      }
      start_obj[2] = WlzAssignObject(start_obj[2], NULL);

      switch( connectivity ){

      case WLZ_8_CONNECTED:
      case WLZ_4_CONNECTED:
	if((dest_obj[1] = WlzDilation(start_obj[1], connectivity,
	                              NULL)) != NULL){
	  dilatobj->domain.p->domains[p] =	
	    WlzAssignDomain(dest_obj[1]->domain, NULL);
	  WlzFreeObj( dest_obj[1] );
	}
	else {
	  dilatobj->domain.p->domains[p].core = NULL;
	}
	break;

      case WLZ_6_CONNECTED:
	dest_obj[0] = WlzAssignObject(start_obj[0], NULL);
	dest_obj[1] = WlzDilation(start_obj[1], WLZ_4_CONNECTED, NULL);
	dest_obj[2] = WlzAssignObject(start_obj[2], NULL);
	tmp_obj = WlzUnionN(3, dest_obj, 0, NULL);
	dilatobj->domain.p->domains[p] = WlzAssignDomain(tmp_obj->domain,
							NULL);
	WlzFreeObj( tmp_obj );
	WlzFreeObj( dest_obj[0] );
	WlzFreeObj( dest_obj[1] );
	WlzFreeObj( dest_obj[2] );
	break;

      case WLZ_18_CONNECTED:
	dest_obj[0] = WlzDilation(start_obj[0], WLZ_4_CONNECTED, NULL);
	dest_obj[1] = WlzDilation(start_obj[1], WLZ_8_CONNECTED, NULL);
	dest_obj[2] = WlzDilation(start_obj[2], WLZ_4_CONNECTED, NULL);
	tmp_obj = WlzUnionN(3, dest_obj, 0, NULL);
	dilatobj->domain.p->domains[p] = WlzAssignDomain(tmp_obj->domain,
							NULL);
	WlzFreeObj( tmp_obj );
	WlzFreeObj( dest_obj[0] );
	WlzFreeObj( dest_obj[1] );
	WlzFreeObj( dest_obj[2] );
	break;

      case WLZ_26_CONNECTED:
	dest_obj[0] = WlzDilation(start_obj[0], WLZ_8_CONNECTED, NULL);
	dest_obj[1] = WlzDilation(start_obj[1], WLZ_8_CONNECTED, NULL);
	dest_obj[2] = WlzDilation(start_obj[2], WLZ_8_CONNECTED, NULL);
	tmp_obj = WlzUnionN(3, dest_obj, 0, NULL);
	dilatobj->domain.p->domains[p] = WlzAssignDomain(tmp_obj->domain,
							NULL);
	WlzFreeObj( tmp_obj );
	WlzFreeObj( dest_obj[0] );
	WlzFreeObj( dest_obj[1] );
	WlzFreeObj( dest_obj[2] );
	break;

      default:
	errNum = WLZ_ERR_CONNECTIVITY_TYPE;
        break;
      }
    }
    WlzFreeObj( start_obj[0] );
    WlzFreeObj( start_obj[1] );
    WlzFreeObj( start_obj[2] );

    WlzStandardPlaneDomain(dilatobj->domain.p, NULL);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return dilatobj;
}
