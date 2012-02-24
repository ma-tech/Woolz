#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzBndJavaUtils_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlzBnd/WlzBndJavaUtils.c
* \author       Nick Burton
* \date         August 2003
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
* \brief	Binding functions to enable Java Woolz to access the
*		Woolz data structures.
* \ingroup	LibWlzBnd
*/
#include <Wlz.h>

/*---------------------------------*/
int WlzBndObjGetType(WlzObject *obj) {
   return obj->type;
}
/*---------------------------------*/
int WlzBndGetIntervalDomainLine1(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_2D_DOMAINOBJ) return -2;
   return obj->domain.i->line1;
}
int WlzBndGetIntervalDomainLastLn(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_2D_DOMAINOBJ) return -2;
   return obj->domain.i->lastln;
}
int WlzBndGetIntervalDomainKol1(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_2D_DOMAINOBJ) return -2;
   return obj->domain.i->kol1;
}
int WlzBndGetIntervalDomainLastKl(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_2D_DOMAINOBJ) return -2;
   return obj->domain.i->lastkl;
}
/*---------------------------------*/
int WlzBndGetPlaneDomainPlane1(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_3D_DOMAINOBJ) return -2;
   return obj->domain.p->plane1;
}
int WlzBndGetPlaneDomainLastPl(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_3D_DOMAINOBJ) return -2;
   return obj->domain.p->lastpl;
}
int WlzBndGetPlaneDomainLine1(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_3D_DOMAINOBJ) return -2;
   return obj->domain.p->line1;
}
int WlzBndGetPlaneDomainLastLn(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_3D_DOMAINOBJ) return -2;
   return obj->domain.p->lastln;
}
int WlzBndGetPlaneDomainKol1(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_3D_DOMAINOBJ) return -2;
   return obj->domain.p->kol1;
}
int WlzBndGetPlaneDomainLastKl(WlzObject *obj) {
   if(obj == NULL) return -1;
   if(obj->type != WLZ_3D_DOMAINOBJ) return -2;
   return obj->domain.p->lastkl;
}
/*---------------------------------*/
WlzObject *WlzBndSetGreyValues(WlzObject           *domainObj,
			   WlzObject           *valuesObj,
			   WlzErrorNum         *dstErr) {

   WlzObject *retObj = NULL;
   WlzErrorNum  errNum = WLZ_ERR_NONE;

   if(domainObj->type != WLZ_2D_DOMAINOBJ) return (NULL);
   if(valuesObj->type != WLZ_2D_DOMAINOBJ) return (NULL);

   retObj = WlzMakeMain(domainObj->type,
		  domainObj->domain,
		  valuesObj->values,
		  NULL,
		  NULL,
		  &errNum);

  *dstErr = errNum;
  return (retObj);

}
/*---------------------------------*/
WlzAffineTransform *WlzBndGetAffineTransform(
                           WlzThreeDViewStruct *viewStr,
			   WlzErrorNum         *dstErr) {

   WlzAffineTransform *retObj = NULL;
   WlzErrorNum  errNum = WLZ_ERR_NONE;

   if(viewStr == NULL) return (NULL);

   retObj = viewStr->trans;

   *dstErr = errNum;
   return (retObj);

}
/*---------------------------------*/
void WlzBndGetTransformMatrix(
               WlzIVertex2 *dstSizeArrayDat,
	       double  ***dstArrayDat,
	       WlzAffineTransform *trans,
	       WlzErrorNum *dstErr) {

   WlzErrorNum   errNum = WLZ_ERR_NONE;
   WlzIVertex2 size;

   size.vtX = 4;
   size.vtY = 4;

   if(dstArrayDat == NULL) {
     errNum = WLZ_ERR_PARAM_NULL;
   } else {
     *dstSizeArrayDat = size;
     *dstArrayDat = trans->mat;
   }

   *dstErr = errNum;
}
/*---------------------------------*/
