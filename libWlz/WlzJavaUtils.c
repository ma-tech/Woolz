#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzJavaUtils.c
* \author       Nick Burton
* \date         September 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
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
int WlzObjGetType(WlzObject *obj) {
   return obj->type;
}
/*---------------------------------*/
int WlzGetIntervalDomainLine1(WlzObject *obj) {
   return obj->domain.i->line1;
}
int WlzGetIntervalDomainLastLn(WlzObject *obj) {
   return obj->domain.i->lastln;
}
int WlzGetIntervalDomainKol1(WlzObject *obj) {
   return obj->domain.i->kol1;
}
int WlzGetIntervalDomainLastKl(WlzObject *obj) {
   return obj->domain.i->lastkl;
}
/*---------------------------------*/
int WlzGetPlaneDomainPlane1(WlzObject *obj) {
   return obj->domain.p->plane1;
}
int WlzGetPlaneDomainLastPl(WlzObject *obj) {
   return obj->domain.p->lastpl;
}
int WlzGetPlaneDomainLine1(WlzObject *obj) {
   return obj->domain.p->line1;
}
int WlzGetPlaneDomainLastLn(WlzObject *obj) {
   return obj->domain.p->lastln;
}
int WlzGetPlaneDomainKol1(WlzObject *obj) {
   return obj->domain.p->kol1;
}
int WlzGetPlaneDomainLastKl(WlzObject *obj) {
   return obj->domain.p->lastkl;
}
/*---------------------------------*/
