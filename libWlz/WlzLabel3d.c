#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzLabel3d.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Segments a 3D Woolz object. Currently 3D connectivity
*		is defined by intersection (ie 6-connected).
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>
#include <Wlz.h>


static WlzObject	***objarray;
static int 		*narray;

static WlzObject *find3Dobj(
  WlzObject		*obj,
  int			p,
  WlzPlaneDomain	*pladm);

/************************************************************************
*   Function   : WlzLabel3d						*
*   Date       : Mon Nov 25 15:00:01 1996				*
*************************************************************************
*   Synopsis   :Segment a 3D object					*
*   Returns    :WlzErrorNum: WLZ_ERR_NONE on success, otherwise:	*
*		WLZ_ERR_PLANEDOMAIN_TYPE, WLZ_ERR_MEM_ALLOC.		*
*   Parameters :WlzObject	*obj: object to be segmented		*
*		int	*mm: number of objects return			*
*		WlzObject	**objlist: object list array for object	*
*			pointer return					*
*		int	nobj: maximum number of objects in array	*
*		int	ignln: ignore objects with num_lines <= ignln	*
*		WlzConnectType connect: connectivity type 4 or 8	*
*   Global refs:None.							*
************************************************************************/

WlzErrorNum WlzLabel3d(WlzObject	*obj,
		       int		*numobj,
		       WlzObject	**objlist,
		       int		nobj,
		       int		ignlns,
		       WlzConnectType	connect)
{
  /* local variables */
  WlzObject		tobj;
  WlzPlaneDomain	*planedm;
  WlzDomain		*domains;
  WlzValues		*values;
  int 			i, j, nplanes;
  WlzErrorNum		wlzerrno=WLZ_ERR_NONE;

  /* only need to check the planedomain type because object pointer
     checks have been done by WlzLabel. This procedure must only be
     called via WlzLabel unless the corresponding checks have been
     made. Do not enter this procedure in WlzProto.h */
  if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
    *numobj = 0;
    return WLZ_ERR_PLANEDOMAIN_TYPE;
  }

  planedm = obj->domain.p;
  domains = planedm->domains;
  if( obj->values.core != NULL ){
    values = obj->values.vox->values;
  }
  else{
    values = NULL;
  }
  nplanes = planedm->lastpl - planedm->plane1 + 1;

  /* allocate space for working arrays */
  if( (objarray = (WlzObject ***)
       AlcMalloc(sizeof(WlzObject **) * nplanes)) == NULL ){
    *numobj = 0;
    return WLZ_ERR_MEM_ALLOC;
  }
  if( (narray = (int *) AlcMalloc(sizeof(int) * nplanes)) == NULL ){
    AlcFree((void *) objarray);
    *numobj = 0;
    return WLZ_ERR_MEM_ALLOC;
  }

  /* find object list at each plane */
  tobj.type = WLZ_2D_DOMAINOBJ;
  tobj.linkcount = 0;
  tobj.plist = NULL;
  tobj.assoc = NULL;
  for(i=0; i < nplanes; i++){
    tobj.domain = domains[i];
    tobj.values.core = NULL;
    if( values ){
      tobj.values = values[i];
    }

    if( (wlzerrno = WlzLabel(&tobj, &narray[i], &objarray[i], nobj,
			     ignlns, connect)) != WLZ_ERR_NONE ){
      /* clean up in here */
      *numobj = 0;
	return wlzerrno;
    }

  }

  /* find sets of objects overlapping in 3D */
  *numobj = 0;
  for(i=0; i < nplanes; i++){
    while((narray[i] != 0 ) && (*numobj < nobj)){
      *numobj += 1;
      objlist[*numobj - 1] = find3Dobj(objarray[i][0], i, planedm);
    }
  }

  /* free allocated space */
  AlcFree((void *) narray);
  for(i=0; i < nplanes; i++){
    if( objarray[i] != NULL ){
      AlcFree((void *) objarray[i]);
    }
  }
  AlcFree((void *) objarray);

  return wlzerrno;
}


static WlzObject *find3Dobj(
  WlzObject		*obj,
  int			p,
  WlzPlaneDomain	*pladm)
{
  /* local variables */
  WlzObject 		*new_3D, *temp_3D, *current_3D;
  WlzObject		*iobj, **objlist, *objn[2];
  WlzPlaneDomain	*pdom;
  WlzVoxelValues	*voxtab;
  WlzDomain		domain;
  WlzValues		values;
  int			i, j, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* remove current object from the array */
  objlist = objarray[p];
  for(i=0, j=0; i < narray[p]; i++, j++){
    if( objlist[i] == obj ){
      j++;
    }
    objlist[i] = objlist[j];
  }

  narray[p]--;

  /* make a 3D object with one plane */
  pdom = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,p,p,
			    obj->domain.i->line1,obj->domain.i->lastln,
			    obj->domain.i->kol1,obj->domain.i->lastkl, 
			    &errNum);
  pdom->domains[0] = WlzAssignDomain( obj->domain, NULL );
  pdom->voxel_size[0] = pladm->voxel_size[0];
  pdom->voxel_size[1] = pladm->voxel_size[1];
  pdom->voxel_size[2] = pladm->voxel_size[2];
  if( obj->values.core ){
    voxtab = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
				 p, p,
				 WlzGetBackground(obj, NULL),
				 NULL, &errNum);
    voxtab->values[0] = WlzAssignValues( obj->values, NULL );
  }
  else {
    voxtab = NULL;
  }
  
  domain.p = pdom;
  values.vox = voxtab;
  current_3D = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			   NULL, NULL, &errNum);

  /* check for overlaps in the plane below */
  if( p != 0 ){
    for(i=0; i < narray[p-1]; i++){
      iobj = WlzIntersect2(obj, objarray[p-1][i], &errNum);
      if( (iobj != NULL) && (WlzIntervalCount(iobj->domain.i, NULL) != 0) ){
	new_3D = find3Dobj(objarray[p-1][i], p-1, pladm);
	objn[0] = current_3D;
	objn[1] = new_3D;
	temp_3D = WlzUnionN(2, &objn[0], 1, &errNum);
	WlzFreeObj(current_3D);
	WlzFreeObj(new_3D);
	current_3D = temp_3D;
      }
      if( iobj != NULL ){
	WlzFreeObj(iobj);
      }
    }
  }

  /* check for overlaps in the plane above */
  nplanes = pladm->lastpl - pladm->plane1 + 1;
  if( p < (nplanes-1) ){
    for(i=0; i < narray[p+1]; i++){
      iobj = WlzIntersect2(obj, objarray[p+1][i], &errNum);
      if( (iobj != NULL) && (WlzIntervalCount(iobj->domain.i, NULL) != 0) ){
	new_3D = find3Dobj(objarray[p+1][i], p+1, pladm);
	objn[0] = current_3D;
	objn[1] = new_3D;
	temp_3D = WlzUnionN(2, objn, 1, &errNum);
	WlzFreeObj(current_3D);
	WlzFreeObj(new_3D);
	current_3D = temp_3D;
      }
      WlzFreeObj(iobj);
    }
  }

  /* return current 3D object, also free the given object! - not a good
     design */
  WlzFreeObj( obj );

  return(current_3D);
}
