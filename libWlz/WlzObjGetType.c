#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzObjGetType.c
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
* \brief	for java woolz; to allow java progs to get type
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>

int WlzObjGetType(WlzObject *obj) {

   return obj->type;
}
