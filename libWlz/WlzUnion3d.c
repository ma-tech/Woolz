#pragma ident "MRC HGU $Id$"
/***********************************************************************
* Project:      Woolz
* Title:        WlzUnion3d.c
* Date:         March 1999
* Author:       Richard Baldock
* Copyright:	1999 Medical Research Council, UK.
*		All rights reserved.
* Address:	MRC Human Genetics Unit,
*		Western General Hospital,
*		Edinburgh, EH4 2XU, UK.
* Purpose:      Computes the set union of 3D Woolz objects.
* $Revision$
* Maintenance:	Log changes below, with most recent at top of list.
************************************************************************/
#include <stdlib.h>

#include <Wlz.h>

/************************************************************************
*   Function   : WlzUnion3d						*
*   Date       : Mon Oct 21 17:06:15 1996				*
*************************************************************************
*   Synopsis   :Private routine for use by WlzUnionN			*
*   Returns    :							*
*   Parameters :							*
*   Global refs:							*
************************************************************************/

WlzObject *WlzUnion3d(int	n,
		      WlzObject **objs,
		      int	uvt,
		      WlzErrorNum *dstErr)
{
  /* local variables */
  WlzObject 		**objlist, *newobj;
  WlzPlaneDomain 	*pdom, *newpdom;
  WlzVoxelValues 	*voxtab, *newvoxtab;
  WlzDomain 		*domains, domain;
  WlzValues	 	*values, vals;
  int 			i, p, np, min_plane, max_plane;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

    /* all objects have been checked by WlzUnionN therefore do not need
       checking here. This routine should not be used except through
       WlzUnionN and must not be included in WlzProto.h */

    /* we do need to check the planedomain types however - this
       used to strip out the wrong types, now it will return an error */
  for (i=0; i<n ; i++ ){
    if ( objs[i]->domain.p->type != WLZ_PLANEDOMAIN_DOMAIN ){
      errNum = WLZ_ERR_DOMAIN_TYPE;
      break;
    }
  }
  newobj = NULL;

  if( (errNum == WLZ_ERR_NONE) && (n == 1) ){
    return WlzMakeMain(objs[0]->type, objs[0]->domain,
		       objs[0]->values, NULL, NULL, dstErr);
  }

  /* find the minimum and maximum plane */
  if( errNum == WLZ_ERR_NONE ){
    pdom = (objs[0])->domain.p;
    min_plane = pdom->plane1;
    max_plane = pdom->lastpl;
    for(i=1; i < n; i++)
    {
      pdom = objs[i]->domain.p;

      if( pdom->plane1 < min_plane )
      {
	min_plane = pdom->plane1;
      }

      if( pdom->lastpl > max_plane )
      {
	max_plane = pdom->lastpl;
      }

      if( objs[i]->values.core == NULL )
      {
	uvt = 0;
      }
    }
  }

  /* allocate space for a working object array */
  objlist = NULL;
  if( errNum == WLZ_ERR_NONE ){
    if( (objlist = (WlzObject **) AlcMalloc(sizeof(WlzObject *) * n))
       == NULL){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      domain.core = NULL;
      vals.core = NULL;
      for(i=0; i < n; i++){
	objlist[i] = WlzMakeMain(WLZ_2D_DOMAINOBJ, domain, vals,
				 NULL, NULL, &errNum);
      }
    }
  }

  /* make a new planedomain and valuetable if required */
  if((errNum == WLZ_ERR_NONE) &&
     (newpdom = WlzMakePlaneDomain(pdom->type, min_plane, max_plane,
				   0, 0, 0, 0, &errNum)) ){
    domains = newpdom->domains;
    newvoxtab = NULL;
    if( uvt ){
      if( newvoxtab = WlzMakeVoxelValueTb((*objs)->values.vox->type,
					  min_plane, max_plane,
					  (*objs)->values.vox->bckgrnd,
					  NULL, &errNum) ){
	values = newvoxtab->values;
      }
    }
  }

  /* find union at each plane */
  if( errNum == WLZ_ERR_NONE ){
    for(p=min_plane; p <= max_plane; p++){
      np = 0;
      for(i=0; i < n; i++){
	pdom = objs[i]->domain.p;
	if( pdom->plane1 > p || pdom->lastpl < p )
	{
	  continue;
	}

	if( (pdom->domains)[p - pdom->plane1].i == NULL )
	{
	  continue;
	}

	objlist[np]->domain.i = (pdom->domains)[p - pdom->plane1].i;
	if( uvt ){
	  voxtab = objs[i]->values.vox;
	  objlist[np]->values.v = (voxtab->values)[p - voxtab->plane1].v;
	}
	np++;
      }

      newobj = NULL;
      if( np ){
	newobj = WlzUnionN(np, objlist, uvt, &errNum);
      }

      /* when np is 1, WlzUnionN does not return a copy. */
      if( newobj != NULL ){
	if (1 == np)
	  domains[p - min_plane] = WlzCopyDomain(newobj->type, newobj->domain, NULL);
	else
	  domains[p - min_plane] = WlzAssignDomain(newobj->domain, NULL);
	if( uvt ){
	  if (1 == np)
	    values[p - min_plane] = WlzCopyValues(newobj->type, newobj->values, newobj->domain, NULL);
	  else
	    values[p - min_plane] = WlzAssignValues(newobj->values, NULL);
	}
	WlzFreeObj(newobj);
      } else {
	domains[p - min_plane].i = NULL;
	if( uvt ){
	  values[p - min_plane].v = NULL;
	}
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    /* copy voxel sizes from first object to new object */
    pdom = (WlzPlaneDomain *) (objs[0])->domain.p;
    newpdom->voxel_size[0] = pdom->voxel_size[0];
    newpdom->voxel_size[1] = pdom->voxel_size[1];
    newpdom->voxel_size[2] = pdom->voxel_size[2];

    /* make a new object */
    WlzStandardPlaneDomain(newpdom, newvoxtab);
    domain.p = newpdom;
    vals.vox = newvoxtab;
    newobj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, vals,
			 NULL, NULL, &errNum);

  /* free allocated memory and return */
    for(i=0; i < n; i++){
      AlcFree(objlist[i]);
    }
    AlcFree( objlist );
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( newobj );
}
