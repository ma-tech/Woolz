#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzIntersect3d.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Calculates the intersection between a list of 3D
*		Woolz objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>

/************************************************************************
*   Function   : WlzIntersect3d						*
*   Date       : Mon Oct 21 19:11:53 1996				*
*************************************************************************
*   Synopsis   :Calculate the intersection between a list of 3D objects.*
*   Returns    :WlzObject *: the intersection, NULL on error. Nil	*
*		intersection is returned as WLZ_EMPTY_OBJ		*
*   Parameters :WlzObject	**objs: array of object pointers.	*
*		int 		n: number of objects.			*
*		int		uvt: copy grey-values flag.		*
*   Global refs:None.							*
************************************************************************/

WlzObject *WlzIntersect3d(WlzObject	**objs,
			  int 		n,
			  int		uvt,
			  WlzErrorNum   *wlzErr)
{
  /* local variables */
  WlzObject 		**objlist, *newObj;
  WlzPlaneDomain 	*pdom, *newpdom;
  WlzVoxelValues	*voxtab, *newvoxtab;
  WlzDomain 		*domains, domain;
  WlzValues 		*values, vals;
  WlzPixelV		bgd;
  int 			i, p, np, min_plane, max_plane, emptyFlag;
  WlzErrorNum		errNum = WLZ_ERR_NONE;

  /* all objects have been checked by WlzIntersectN therefore do not need
     checking here. This routine should not be used except through
     WlzIntersectN and must not be included in WlzProto.h */

  /* check all objects are non-empty and have the same type
     Note an empty object is not an error */
  for (i=0; i<n; i++){
    if( objs[i]->type != objs[0]->type ){
      if( objs[i]->type == WLZ_EMPTY_OBJ ){
	return WlzMakeEmpty(wlzErr);
      }
      else {
	newObj = NULL;
	errNum = WLZ_ERR_OBJECT_TYPE;
	if(wlzErr) {
	  *wlzErr = errNum;
	}
	return(newObj);
      }
    }

    /* check for size */
    if( WlzIsEmpty(objs[i], &errNum) ){
      return WlzMakeEmpty(wlzErr);
    }
    else {
      if( errNum != WLZ_ERR_NONE ){
	if(wlzErr) {
	  *wlzErr = errNum;
	}
	return NULL;
      }
    }
  }

  /* we do need to check the planedomain types */
  for(i=0; i < n; i++){
    if( objs[i]->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
      if(wlzErr) {
        *wlzErr = WLZ_ERR_DOMAIN_TYPE;
      }
      return NULL;
    }
  }

  /* check number */
  if (n == 1){
    newObj = WlzMakeMain(objs[0]->type, objs[0]->domain, objs[0]->values,
			 NULL, NULL, &errNum);
    if(wlzErr) {
      *wlzErr = errNum;
    }
    return(newObj);
  }

  /* find the minimum and maximum plane */
  pdom = objs[0]->domain.p;
  min_plane = pdom->plane1;
  max_plane = pdom->lastpl;
  for(i=1; i < n; i++){
    pdom = objs[i]->domain.p;

    if( pdom->plane1 > min_plane ){
      min_plane = pdom->plane1;
    }

    if( pdom->lastpl < max_plane ){
      max_plane = pdom->lastpl;
    }

    if( objs[i]->values.v == NULL ){
      uvt = 0;
    }
  }
  if( min_plane > max_plane ){
    return WlzMakeEmpty(wlzErr);
  }

  /* allocate space for a working object array */
  if( (objlist = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * n)) == NULL){
    if(wlzErr) {
      *wlzErr = WLZ_ERR_MEM_ALLOC;
    }
    return NULL;
  }
  domain.core = NULL;
  vals.core = NULL;
  for(i=0; (i < n) && (errNum == WLZ_ERR_NONE); i++){
    objlist[i] = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, vals, NULL, NULL,
    			     &errNum);
  }
  if(errNum != WLZ_ERR_NONE) {
    if(wlzErr) {
      *wlzErr = errNum;
    }
    return NULL;
  }

  /* make a new planedomain and valuetable if required */
  newpdom = WlzMakePlaneDomain(pdom->type, min_plane, max_plane, 0, 0, 0, 0,
  			       &errNum);
  if(errNum == WLZ_ERR_NONE) {
    domains = newpdom->domains;
    newvoxtab = NULL;
    if( uvt ){
      bgd.type = WLZ_GREY_INT;
      bgd.v.inv = 0;
      newvoxtab = WlzMakeVoxelValueTb((*objs)->values.vox->type, min_plane,
				      max_plane, bgd, NULL, &errNum);
      values = newvoxtab->values;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    if(wlzErr) {
      *wlzErr = errNum;
    }
    return NULL;
  }

  /* find intersection at each plane */
  for(p=min_plane; (p <= max_plane) && (errNum == WLZ_ERR_NONE); p++){
    np = 0;
    for(i=0; i < n; i++){
      pdom = objs[i]->domain.p;

      objlist[np]->domain.i = (pdom->domains)[p - pdom->plane1].i;
      if( uvt ){
	voxtab = objs[i]->values.vox;
	objlist[np]->values.v = (voxtab->values)[p - voxtab->plane1].v;
      }
      if( objlist[np]->domain.i ){
	np++;
      }
    }

    newObj = NULL;
    if( np == n ){
      newObj = WlzIntersectN(np, objlist, uvt, &errNum);
    }
    else {
      newObj = WlzMakeEmpty(NULL);
    }
    if(newObj) {
      emptyFlag = (newObj->type == WLZ_EMPTY_OBJ);
    }
    else {
      if(wlzErr) {
	*wlzErr = errNum;
      }
      return NULL;
    }

    if(emptyFlag){
      domains[p - min_plane].core = NULL;
      if( uvt ){
	values[p - min_plane].core = NULL;
      }
    }
    else {
      domains[p - min_plane] = WlzAssignDomain(newObj->domain, &errNum);




      if(uvt && (errNum == WLZ_ERR_NONE)){
	values[p - min_plane] = WlzAssignValues(newObj->values, &errNum);
      }
    }
    WlzFreeObj(newObj);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* copy voxel sizes from first object to new object */
    pdom = objs[0]->domain.p;
    newpdom->voxel_size[0] = pdom->voxel_size[0];
    newpdom->voxel_size[1] = pdom->voxel_size[1];
    newpdom->voxel_size[2] = pdom->voxel_size[2];

      /* make a new object */
    (void )WlzStandardPlaneDomain(newpdom, newvoxtab);
    domain.p = newpdom;	
    vals.vox = newvoxtab;
    newObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, vals, NULL, NULL, &errNum);
  }
  else {
    newObj = NULL;
  }

  /* free allocated memory and return */
  for(i=0; i < n; i++){
    AlcFree(objlist[i]);
  }
  AlcFree( objlist );
  if(wlzErr) {
    *wlzErr = errNum;
  }

  return( newObj );
}
