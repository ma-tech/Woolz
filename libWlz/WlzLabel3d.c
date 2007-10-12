#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzLabel3d_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         libWlz/WlzLabel3d.c
* \author       Richard Baldock
* \date         August 2003
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
* \brief	Segments a 3D object, called from WlzLabel().
* \ingroup	WlzBinaryOps
* \todo         -
* \bug          None known.
*/

#include <stdlib.h>
#include <Wlz.h>


static WlzObject	***objarray;
static int 		*narray;

static WlzObject *find3Dobj(
  WlzObject		*obj,
  int			p,
  WlzPlaneDomain	*pladm,
  WlzErrorNum		*dstErr);

static WlzErrorNum removeObj(
  WlzObject	*obj,
  int		p);

/* function:     WlzLabel3d    */
/*! 
* \ingroup      WlzBinaryOps
* \brief        Segment a 3D object. This is a private routine for
 WlzLabel() and should not be called directly.
*
* \return       Error number.
* \param    obj	Input object to be segmented.
* \param    numobj	Number of objects found.
* \param    objlist	Array of object pointers.
* \param    nobj	Maximum number of objects in array.
* \param    ignlns	ignore objects with num. line <= ignlns.
* \param    connect	Connectivity of segmented objects.
* \par      Source:
*                WlzLabel3d.c
*/
WlzErrorNum WlzLabel3d(
  WlzObject	*obj,
  int		*numobj,
  WlzObject	**objlist,
  int		nobj,
  int		ignlns,
  WlzConnectType	connect)
{
  /* local variables */
  WlzObject		*tobj;
  WlzPlaneDomain	*planedm;
  WlzDomain		*domains;
  WlzValues		*values;
  int 			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* only need to check the planedomain type because object pointer
     checks have been done by WlzLabel. This procedure must only be
     called via WlzLabel unless the corresponding checks have been
     made. Do not enter this procedure in WlzProto.h */
  if( obj->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
    *numobj = 0;
    return WLZ_ERR_PLANEDOMAIN_TYPE;
  }

  /* check the connectivity */
  switch( connect ){
  case WLZ_4_CONNECTED:
  case WLZ_6_CONNECTED:
    connect = WLZ_4_CONNECTED;
    break;

  case WLZ_8_CONNECTED:
  case WLZ_18_CONNECTED:
  case WLZ_26_CONNECTED:
  default:
    connect = WLZ_8_CONNECTED;
    break;
  }

  /* set up local variables */
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
  if( errNum == WLZ_ERR_NONE ){
    for(i=0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++){
      WlzDomain	domain;
      WlzValues	value;

      domain = domains[i];
      value.core = NULL;
      if( values ){
	value = values[i];
      }
      if((domain.core != NULL) && (value.core != NULL) &&
         ((tobj = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, value,
			      NULL, NULL, &errNum)) != NULL)){
	tobj = WlzAssignObject(tobj, NULL);
	errNum = WlzLabel(tobj, &narray[i], &objarray[i], nobj,
			    ignlns, connect);
	WlzFreeObj(tobj);
      }
      else
      {
        narray[i] = 0;
      }
    }
  }

  /* find sets of objects overlapping in 3D */
  if( errNum == WLZ_ERR_NONE ){
    *numobj = 0;
    for(i=0; (i < nplanes) && (errNum == WLZ_ERR_NONE); i++){
      while((narray[i] != 0 ) && (*numobj < nobj)){
	*numobj += 1;
	tobj = objarray[i][0];
	errNum = removeObj(tobj, i);
	objlist[*numobj - 1] = find3Dobj(tobj, i, planedm, &errNum);
	WlzFreeObj(tobj);
      }
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

  return errNum;
}

static WlzErrorNum removeObj(
  WlzObject	*obj,
  int		p)
{
  WlzObject	**objlist;
  int		i, j;

  /* remove object from the global array */
  objlist = objarray[p];
  for(i=0, j=0; i < narray[p]; i++, j++){
    if( objlist[i] == obj ){
      narray[p]--;
      j++;
    }
    if( i < narray[p] ){
      objlist[i] = objlist[j];
    }
  }

  if( i == j ){
    return WLZ_ERR_UNSPECIFIED;
  }
  return WLZ_ERR_NONE;
}



static WlzObject *find3Dobj(
  WlzObject		*obj,
  int			p,
  WlzPlaneDomain	*pladm,
  WlzErrorNum		*dstErr)
{
  /* local variables */
  WlzObject 		*new_3D, *temp_3D, *current_3D, *tobj;
  WlzObject		*objn[2];
  WlzPlaneDomain	*pdom;
  WlzVoxelValues	*voxtab;
  WlzDomain		domain;
  WlzValues		values;
  int			i, nplanes;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  /* make a 3D object with one plane */
  pdom = WlzMakePlaneDomain(WLZ_PLANEDOMAIN_DOMAIN,
			    pladm->plane1 + p, pladm->plane1 + p,
			    obj->domain.i->line1,obj->domain.i->lastln,
			    obj->domain.i->kol1,obj->domain.i->lastkl, 
			    &errNum);
  pdom->domains[0] = WlzAssignDomain( obj->domain, NULL );
  pdom->voxel_size[0] = pladm->voxel_size[0];
  pdom->voxel_size[1] = pladm->voxel_size[1];
  pdom->voxel_size[2] = pladm->voxel_size[2];
  if( obj->values.core ){
    voxtab = WlzMakeVoxelValueTb(WLZ_VOXELVALUETABLE_GREY,
				 pladm->plane1 + p, pladm->plane1 + p,
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
      if( WlzHasIntersection(obj, objarray[p-1][i], &errNum) ){
	tobj = objarray[p-1][i];
	removeObj(tobj, p-1);
	new_3D = find3Dobj(tobj, p-1, pladm, &errNum);
	WlzFreeObj(tobj);
	objn[0] = current_3D;
	objn[1] = new_3D;
	temp_3D = WlzUnionN(2, &objn[0], 1, &errNum);
	WlzFreeObj(current_3D);
	WlzFreeObj(new_3D);
	current_3D = temp_3D;
      }
    }
  }

  /* check for overlaps in the plane above */
  nplanes = pladm->lastpl - pladm->plane1 + 1;
  if( p < (nplanes-1) ){
    for(i=0; i < narray[p+1]; i++){
      if( WlzHasIntersection(obj, objarray[p+1][i], &errNum) ){
	tobj = objarray[p+1][i];
	removeObj(tobj, p+1);
	new_3D = find3Dobj(tobj, p+1, pladm, &errNum);
	WlzFreeObj(tobj);
	objn[0] = current_3D;
	objn[1] = new_3D;
	temp_3D = WlzUnionN(2, objn, 1, &errNum);
	WlzFreeObj(current_3D);
	WlzFreeObj(new_3D);
	current_3D = temp_3D;
      }
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return(current_3D);
}
