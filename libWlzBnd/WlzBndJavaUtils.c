#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBndJavaUtils_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzBndJavaUtils.c
* \author       Nick Burton
* \date         August 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	functions to enable java woolz to access
*		Wlz data structures
* \todo         -
* \bug          None known.
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
