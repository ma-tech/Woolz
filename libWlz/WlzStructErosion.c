#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzStructErosion_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzStructErosion.c
* \author       Richard Baldock
* \date         March 1999
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
* \brief	Performs erosion using a structuring element.
* \ingroup	WlzMorphologyOps
*/

#include <stdlib.h>
#include <Wlz.h>

static int 			intersecitvs(
				  WlzIntervalLine *itva,
				  WlzIntervalLine *itvb,
				  int n,
				  WlzInterval *jp,
				  WlzInterval *aa,
				  WlzInterval *bb,
				  WlzInterval *cc);

static int 			line_struct_shrink(
				  WlzIntervalLine *itvl,
				  WlzInterval *jntl,
				  WlzInterval *ll);

static int 			arr_arr_intersec(
				  WlzInterval *aa,
				  int m,
				  WlzInterval *bb,
				  int n,
				  WlzInterval *cc);

static int 			intl_to_intl(
				  WlzInterval *intl,
				  WlzInterval *jntl,
				  int n);

static WlzObject 		*WlzStructErosion3d(
				  WlzObject *obj,
				  WlzObject *structElm,
				  WlzErrorNum *dstErr);

/*!
* \return	New object or NULL on error.
* \ingroup 	WlzMorphologyOps
* \brief	Performs erosion using a structuring element.
* \param	obj			Given object to be eroded.
* \param	structElm		Structuring element.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzStructErosion(
  WlzObject 	*obj,
  WlzObject	*structElm,
  WlzErrorNum	*dstErr)
{
  WlzObject		*rtnObj=NULL;
  WlzDomain		bDom,
  			sDom,
  			domain;
  WlzValues		values;
  WlzInterval		*jp, *jpw;
  WlzIntervalLine	*bitv, *sitv;
  WlzInterval 		*aa = NULL;
  WlzInterval 		*bb = NULL;
  WlzInterval 		*cc = NULL;
  int			i, j, k, m;
  int			maxItvLn;
  int			line1, kol1, lastln, lastkl;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  bDom.core = NULL;
  sDom.core = NULL;
  /* check the object */
  if( obj == NULL ){
    errNum = WLZ_ERR_OBJECT_NULL;
  }
  else {
    switch( obj->type ){
    case WLZ_2D_DOMAINOBJ:
      /* check the domain */
      if( obj->domain.core == NULL ){
	errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else if( obj->domain.i->type == WLZ_EMPTY_DOMAIN ){
	return WlzMakeEmpty(dstErr);
      }
      break;

    case WLZ_3D_DOMAINOBJ:
      return WlzStructErosion3d(obj, structElm, dstErr);

    case WLZ_TRANS_OBJ:
      if((values.obj = WlzStructErosion(obj->values.obj, structElm,
					&errNum)) != NULL){
	return WlzMakeMain(WLZ_TRANS_OBJ, obj->domain, values,
			   NULL, NULL, dstErr);
      }
      break;

    case WLZ_EMPTY_OBJ:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_OBJECT_TYPE;
      break;
    }
  }

  /* check the structuring element */
  if( errNum == WLZ_ERR_NONE ){
    if( structElm == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      switch( structElm->type ){
      case WLZ_2D_DOMAINOBJ:
	/* check the domain */
	if( structElm->domain.core == NULL ){
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if( structElm->domain.i->type == WLZ_EMPTY_DOMAIN ){
	  WlzMakeMain(obj->type, obj->domain, values,
		      NULL, NULL, dstErr);
	}
	break;

      case WLZ_TRANS_OBJ:
	return WlzStructErosion(obj, structElm->values.obj, &errNum);

      case WLZ_EMPTY_OBJ:
	values.core = NULL;
	return WlzMakeMain(obj->type, obj->domain, values,
			   NULL, NULL, dstErr);

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
	break;
      }
    }
  }

  /* If we get this far we have 2D object and structuring element of
     domain type and with non-null domains */
  if(errNum == WLZ_ERR_NONE)
  {
    if(obj->domain.i->type == WLZ_INTERVALDOMAIN_INTVL)
    {
      bDom = WlzAssignDomain(obj->domain, NULL);
    }
    else
    {
      bDom.i = WlzNewIDomain(WLZ_INTERVALDOMAIN_INTVL, obj->domain.i,
                             &errNum);
      (void )WlzAssignDomain(bDom, NULL);
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(structElm->domain.i->type == WLZ_INTERVALDOMAIN_INTVL)
    {
      sDom = WlzAssignDomain(structElm->domain, NULL);
    }
    else
    {
      sDom.i = WlzNewIDomain(WLZ_INTERVALDOMAIN_INTVL, structElm->domain.i,
                             &errNum);
      (void )WlzAssignDomain(sDom, NULL);
    }
  }

  if(errNum == WLZ_ERR_NONE)
  {
    /* check if the object will be reduced to an empty domain,
       this assumes that every line of the structuring element
       has at least one interval */
    if(sDom.i->lastln - sDom.i->line1 > bDom.i->lastln - bDom.i->line1){
      return WlzMakeEmpty(dstErr);
    }

    /*
     * space for the intervals of the resultant object.
     */
    m = 2 * WlzIntervalCount(obj->domain.i, NULL);
    if((jp = (WlzInterval *)AlcMalloc(sizeof(WlzInterval) * m)) != NULL){
      jpw = jp;
      /*
       * construct an interval domain approximately and the return object.
       * This assumes that the structuring element domain is "standard" ie
       * at least one interval on the first and last lines and the min left
       * offset = 0 and max right offset = width-1.
       */
      line1 = bDom.i->line1 - sDom.i->line1;
      lastln = bDom.i->lastln - sDom.i->lastln;
      kol1 = bDom.i->kol1 - sDom.i->kol1;
      lastkl = bDom.i->lastkl - sDom.i->lastkl;
      if((domain.i = WlzMakeIntervalDomain(WLZ_INTERVALDOMAIN_INTVL,
					   line1, lastln, kol1, lastkl,
					   &errNum)) != NULL){
	domain.i->freeptr = AlcFreeStackPush(domain.i->freeptr, (void *)jp,
					     NULL);
	values.core = NULL;
	rtnObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, values,
			     NULL, NULL, &errNum);
      }
      else {
	AlcFree((void *) jp);
      }
    }
    else {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /* Make buffers with room for the maximum number of intervals in any line. */
  if( errNum == WLZ_ERR_NONE ){
    maxItvLn = 2 * (WlzIDomMaxItvLn(sDom.i) + WlzIDomMaxItvLn(bDom.i));
    if(maxItvLn < 1){
      errNum = WLZ_ERR_DOMAIN_DATA;
    }
    else if(((aa = (WlzInterval *)
                   AlcMalloc(sizeof(WlzInterval) * maxItvLn)) == NULL) ||
            ((bb = (WlzInterval *)
                   AlcMalloc(sizeof(WlzInterval) * maxItvLn)) == NULL) ||
            ((cc = (WlzInterval *)
                   AlcMalloc(sizeof(WlzInterval) * maxItvLn)) == NULL)) {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  /*
   * scan from first line to last line
   */
  if( errNum == WLZ_ERR_NONE ){

    k = sDom.i->lastln - sDom.i->line1 + 1;
    j = 0;
    bitv = &(bDom.i->intvlines[0]);
    sitv = &(sDom.i->intvlines[0]);
    for(i = line1; i <= lastln; i++){
      j = intersecitvs(bitv, sitv, k, jp, &aa[0], &bb[0], &cc[0]);
      WlzMakeInterval(i, domain.i, j, jp);
      jp += j;
      bitv++;
    }
    /*
     * final adjust - check for empty object and standardise the domain
     */
    if(jp == jpw){
      WlzFreeObj(rtnObj);
      rtnObj = WlzMakeEmpty(&errNum);
    }
    else {
      if( (errNum = WlzStandardIntervalDomain(rtnObj->domain.i)) !=
	 WLZ_ERR_NONE ){
	WlzFreeObj(rtnObj);
	rtnObj = NULL;
      }
    }
  }
  if(bDom.core)
  {
    (void )WlzFreeDomain(bDom);
  }
  if(sDom.core)
  {
    (void )WlzFreeDomain(sDom);
  }
  AlcFree(aa);
  AlcFree(bb);
  AlcFree(cc);
  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

/*!
* \return
* \ingroup 	WlzMorphologyOpsitva
* \brief	
* \param	itvb
* \param	n
* \param	jp
* \param	aa
* \param	bb
* \param	cc
*/
static int intersecitvs(
  WlzIntervalLine 	*itva,
  WlzIntervalLine	*itvb,
  int			n,
  WlzInterval 	*jp,
  WlzInterval	*aa,
  WlzInterval	*bb,
  WlzInterval	*cc)
{
  WlzInterval *intl;
  int i, k, m;
  int j, pp;

  if(n == 0)
    return(0);
  /*
   * The result of intervals erosed by a interval put in aa
   * The result of line's intersection is kept in cc and bb in turn
   * depending on even line or odd line
   */
  pp = 0;
  m = line_struct_shrink(itva, itvb->intvs, bb);
  intl_to_intl(bb, cc, m);
  for(i = 0; i < n; i++){
    intl = itvb->intvs;
    for(j = 0; j < itvb->nintvs; j++){
      k = line_struct_shrink(itva, intl, aa);
      if(pp == 0){
	m = arr_arr_intersec(aa, k, bb, m, cc);
	pp = 1;
      } else {
	m = arr_arr_intersec(aa, k, cc, m, bb);
	pp = 0;
      }
      intl++;
    }
    itva++;
    itvb++;
  }
  if(pp == 0){
    for(i = 0; i < m; i++){
      jp->ileft = bb[i].ileft;
      jp->iright = bb[i].iright;
      jp++;
    }
  } else {
    for(i = 0; i < m; i++){
      jp->ileft = cc[i].ileft;
      jp->iright = cc[i].iright;
      jp++;
    }
  }
  return(m);
}

/*!
* \return
* \ingroup      WlzMorphologyOps
* \brief	An interval line itvl shrinks by a factor (length of jntl) and
* 		translates according to the position of jntl.
* \param	itvl
* \param	jntl
* \param	ll
*/
static int line_struct_shrink(
  WlzIntervalLine	*itvl,
  WlzInterval		*jntl,
  WlzInterval		*ll)
{
  WlzInterval *intl = itvl->intvs;
  WlzInterval *ptr;
  int n = itvl->nintvs; 
  int i;

  ptr = ll;
  i = jntl->iright - jntl->ileft - 1;
  while(n-- > 0){
    if(intl->iright - intl->ileft > i){
      ll->ileft = intl->ileft - jntl->ileft;
      ll->iright = intl->iright - jntl->iright;
      ll++;
    } 
    intl++;
  }
  return((int)(ll - ptr));
}

/*!
* \return	Number of intervals in intersection.
* \ingroup	WlzMorphologyOps
* \brief	Find the intersection intervals between two lines.
* \param	aa			One set of intervals.
* \param	m			Number of intervals in aa.
* \param	bb			Other set of intervals.
* \param	n			Number of intervals in bb.
* \param	cc			Buffer for intervals.
*/
static int arr_arr_intersec(
  WlzInterval	*aa,
  int		m,
  WlzInterval	*bb,
  int		n,
  WlzInterval	*cc)
{
  WlzInterval *tmp = cc;

  while(m > 0 && n > 0){
    if(aa->iright < bb->ileft){
      aa++;
      m--;
    } else if(bb->iright < aa->ileft){
      bb++;
      n--;
    } else {
      cc->ileft = (aa->ileft > bb->ileft) ?
	aa->ileft : bb->ileft;
      if(aa->iright < bb->iright){
	cc->iright = aa->iright;
	aa++;
	m--;
      } else {
	cc->iright = bb->iright;
	bb++;
	n--;
      }
      cc++;
    }
  }
  return((int)(cc-tmp));
}

/*!
* \return	Zero.
* \ingroup 	WlzMorphologyOps
* \brief	Assigns the value of intl to jntl.
* \param	intl			First array of intervals.
* \param	jntl			Second array of intervals.
* \param	n			Number of intervals.
*/
static int intl_to_intl(
  WlzInterval	*intl,
  WlzInterval	*jntl,
  int		n)
{
  while(n-- > 0){
    jntl->ileft = intl->ileft;
    jntl->iright = intl->iright;
    jntl++;
    intl++;
  }

  return( 0 );
}

/*!
* \return	New object or NULL on error.
* \ingroup 	WlzMorphologyOps
* \brief	Performs structur
* \param	obj			Given object to be eroded.
* \param	structElm		Structuring element.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzStructErosion3d(
  WlzObject 	*obj,
  WlzObject	*structElm,
  WlzErrorNum	*dstErr)
{
  WlzObject	*rtnObj = NULL;
  WlzObject	*obj1 = NULL, *obj2 = NULL, *obj3 = NULL;
  WlzObject	**objList = NULL;
  WlzDomain	domain, *domains = NULL, *domains1 = NULL, *domains2 = NULL;
  WlzValues	values;
  int		i, j, p, plane1, lastpl, nStructPlanes;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* the object is definitely 3D but the domain needs checking */
  if( obj->domain.p == NULL ){
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else {
    switch( obj->domain.p->type ){
    case WLZ_PLANEDOMAIN_DOMAIN:
      break;

    case WLZ_EMPTY_DOMAIN:
      return WlzMakeEmpty(dstErr);

    default:
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }

  /* check the structuring element, a 2D element is applied to each
     plane by creating a temporary 3D object */
  if( errNum == WLZ_ERR_NONE ){
    if( structElm == NULL ){
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else {
      switch( structElm->type ){
      case WLZ_2D_DOMAINOBJ:
	/* check the domain */
	if( structElm->domain.i == NULL ){
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else if( structElm->domain.i->type == WLZ_EMPTY_DOMAIN ){
	  values.core = NULL;
	  return WlzMakeMain(obj->type, obj->domain, values,
			     NULL, NULL, dstErr);
	}
	else if((domain.p =
		 WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				   0, 0,
				   structElm->domain.i->line1,
				   structElm->domain.i->lastln,
				   structElm->domain.i->kol1,
				   structElm->domain.i->lastkl,
				   &errNum)) != NULL){
	  domain.p->domains[0] = WlzAssignDomain(structElm->domain, NULL);
	  values.core = NULL;
	  if((obj1 = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
				 NULL, NULL, &errNum)) != NULL){
	    rtnObj = WlzStructErosion3d(obj, obj1, &errNum);
	    WlzFreeObj(obj1);
	  }
	  else {
	    WlzFreePlaneDomain(domain.p);
	  }
	}
	break;
	
      case WLZ_3D_DOMAINOBJ:
	/* check the domain and type */
	if( structElm->domain.p == NULL ){
	  errNum = WLZ_ERR_DOMAIN_NULL;
	}
	else {
	  switch( structElm->domain.p->type ){
	  case WLZ_PLANEDOMAIN_DOMAIN:
	    break;

	  case WLZ_EMPTY_DOMAIN:
	    return WlzMakeMain(obj->type, obj->domain, values,
			       NULL, NULL, dstErr);

	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	  }
	}
	break;

      case WLZ_TRANS_OBJ:
	return WlzStructErosion3d(obj, structElm->values.obj, dstErr);

      case WLZ_EMPTY_OBJ:
	values.core = NULL;
	return WlzMakeMain(obj->type, obj->domain, values,
			   NULL, NULL, dstErr);

      default:
	errNum = WLZ_ERR_OBJECT_TYPE;
      }
    }
  }

  /* this far with no errors and no rtnObj then the object and structuring
     element are 3D and have non-null and non-empty plane domains.
     We assume that the structuring element and object are standardised
     */
  if( (errNum == WLZ_ERR_NONE) && !rtnObj ){
    /* make a new planedomain with planes eroded and shifted by the
       planes of the structuring element, leave the lines and kols 
       to be sorted out by WlzStandardPlaneDomain */
    /* assume at least one pixel on each structuring element plane */
    plane1 = obj->domain.p->plane1 - structElm->domain.p->plane1;
    lastpl = obj->domain.p->lastpl - structElm->domain.p->lastpl;
    if( lastpl < plane1 ){
      rtnObj = WlzMakeEmpty(NULL);
    }
  }

  /* now finally do the erosion */
  if( (errNum == WLZ_ERR_NONE) && !rtnObj ){
    domain.p = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
				  plane1, lastpl,
				  obj->domain.p->line1,
				  obj->domain.p->lastln,
				  obj->domain.p->kol1,
				  obj->domain.p->lastkl,
				  &errNum);
    domain.p->voxel_size[0] = obj->domain.p->voxel_size[0];
    domain.p->voxel_size[1] = obj->domain.p->voxel_size[1];
    domain.p->voxel_size[2] = obj->domain.p->voxel_size[2];
    domains = domain.p->domains;
    domains1 = obj->domain.p->domains;
    domains2 = structElm->domain.p->domains;
    values.core = NULL;

    nStructPlanes = structElm->domain.p->lastpl -
      structElm->domain.p->plane1 + 1;
    objList = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * nStructPlanes);

    for(p=plane1; p <= lastpl; p++, domains++){
      for(i=0; i < nStructPlanes; i++){
	j = p + structElm->domain.p->plane1 + i - obj->domain.p->plane1;
	if( domains1[j].core ){
	  obj1 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains1[j], values,
			     NULL, NULL, NULL);
	}
	else {
	  obj1 = WlzMakeEmpty(NULL);
	}
	obj1 = WlzAssignObject(obj1, NULL);
	if( domains2[i].core ){
	  obj2 = WlzMakeMain(WLZ_2D_DOMAINOBJ, domains2[i], values,
			     NULL, NULL, NULL);
	}
	else {
	  obj2 = WlzMakeEmpty(NULL);
	}
	obj2 = WlzAssignObject(obj2, NULL);
	objList[i] = WlzAssignObject(WlzStructErosion(obj1, obj2, NULL),
				     NULL);
	WlzFreeObj(obj1);
	WlzFreeObj(obj2);
      }
      obj3 = WlzIntersectN(nStructPlanes, objList, 0, &errNum);
      if( (obj3 == NULL) || (obj3->type == WLZ_EMPTY_OBJ) ){
	(*domains).core = NULL;
      }
      else {
	*domains = WlzAssignDomain(obj3->domain, NULL);
      }
      if( obj3 ){
	WlzFreeObj(obj3);
      }
      for(i=0; i < nStructPlanes; i++){
	if( objList[i] ){
	  WlzFreeObj(objList[i]);
	}
      }
    }
    AlcFree((void *) objList);
    WlzStandardPlaneDomain(domain.p, NULL);
    values.core = NULL;
    rtnObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			 NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnObj;
}

