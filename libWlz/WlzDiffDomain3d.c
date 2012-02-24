#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDiffDomain3d_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzDiffDomain3d.c
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
* \brief	Functions for computing the domain difference between
* 		3D objects.
* \ingroup	WlzDomainOps
*/

#include <stdlib.h>
#include <Wlz.h>

/*!
* \return	The difference object, NULL on error.
* \ingroup	WlzDomainOps
* \brief	Calculates the domain difference between two 3D objects
*		This function assumes that the object pointers have been
*		checked to be non-NULL, non-empty, both 3D types and with
*		non-NULL domains. This is the case if the function is accessed
*		via WlzDiffDomain which is recommended.
* \param	obj1			First object.
* \param	obj2			Second object.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject *WlzDiffDomain3d(
  WlzObject *obj1,
  WlzObject *obj2,
  WlzErrorNum	*dstErr)
{
  /* local variables */
  WlzObject		*newobj, *obj, o1, o2;
  WlzPlaneDomain	*newpdom, *oldpdom1, *oldpdom2;
  WlzVoxelValues	*newvdom, *oldvdom;
  WlzDomain		domain;
  WlzValues		values;
  int			i, j, p;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

    /* Don't need to check objects because WlzDiffDomain3d is only accessed
       via WlzDiffDomain. The prototype for WlzDiffDomain3d must not be put
       in WlzProto.h and the procedure must not be called directly.
       We do need to check the planedomain type however */
  newobj = NULL;
  if( obj1->domain.core->type != WLZ_2D_DOMAINOBJ ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  if( (errNum == WLZ_ERR_NONE) &&
     (obj1->domain.core->type != obj2->domain.core->type) ){
    errNum = WLZ_ERR_DOMAIN_TYPE;
  }

  /* make a new object with same number of planes as obj1 */
  if( errNum == WLZ_ERR_NONE ){
    oldpdom1 = (WlzPlaneDomain *) obj1->domain.p;
    newpdom = WlzMakePlaneDomain(oldpdom1->type,
				 oldpdom1->plane1, oldpdom1->lastpl,
				 oldpdom1->line1, oldpdom1->lastln,
				 oldpdom1->kol1, oldpdom1->lastkl,
				 &errNum);
    if( newpdom ){
      for(i=0; i < 3; i++){
	newpdom->voxel_size[i] = oldpdom1->voxel_size[i];
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    if( obj1->values.core != NULL ){
      oldvdom = obj1->values.vox;
      newvdom = WlzMakeVoxelValueTb(oldvdom->type, oldvdom->plane1,
				    oldvdom->lastpl,
				    oldvdom->bckgrnd, NULL, &errNum);
      if( newvdom ){
	for(i=0; i < oldvdom->lastpl - oldvdom->plane1 + 1; i++){
	  newvdom->values[i] = WlzAssignValues(oldvdom->values[i], NULL);
	}
      }
    }
    else {
      newvdom = NULL;
    }
  }
    
  /* find the new domains */
  if( errNum == WLZ_ERR_NONE ){
    oldpdom2 = obj2->domain.p;
    o1.type = o2.type = WLZ_2D_DOMAINOBJ;
    o1.linkcount = o2.linkcount = 0;
    o1.values.core = o2.values.core = NULL;
    o1.plist = o2.plist = NULL;
    o1.assoc = o2.assoc = NULL;

    for(p=newpdom->plane1, i=0; p <= newpdom->lastpl; p++, i++){

      /* test for non-overlapping planes and make a copy */
      if( p < oldpdom2->plane1 || p > oldpdom2->lastpl ){
	newpdom->domains[i] =
	  WlzAssignDomain(WlzCopyDomain(oldpdom1->type, oldpdom1->domains[i],
					&errNum), NULL);
	continue;
      }

      /* test for NULL domain in obj1
	 this will become an error - replaced by WLZ_EMPTY_DOMAIN */
      if( oldpdom1->domains[i].core == NULL ){
	newpdom->domains[i].core = NULL;
	continue;
      }

      /* get corresponding planes */
      j = p - oldpdom2->plane1;
      o1.domain = oldpdom1->domains[i];
      o2.domain = oldpdom2->domains[j];

      /* test for NULL domain in obj2
	 this will become an error - replaced by WLZ_EMPTY_DOMAIN */
      if( oldpdom2->domains[j].core == NULL ){
	newpdom->domains[i] = 
	  WlzAssignDomain(WlzCopyDomain(oldpdom1->type, oldpdom1->domains[i],
					&errNum), NULL);
	continue;
      }

      /* find difference domain */
      if((obj = WlzDiffDomain( &o1, &o2, &errNum)) != NULL){
	newpdom->domains[i] = WlzAssignDomain(obj->domain, NULL);
	WlzFreeObj( obj );
      }
      else {
	break;
      }
    }
  }

  if( errNum == WLZ_ERR_NONE ){
    WlzStandardPlaneDomain( newpdom, newvdom );
    domain.p = newpdom;
    values.vox = newvdom;
    newobj = WlzMakeMain(WLZ_3D_DOMAINOBJ, domain, values,
			 NULL, NULL, &errNum);
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return newobj;
}
