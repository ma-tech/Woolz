#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzPolyUtils.c
* \author       Bill Hill
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
* \brief	Functions for manipulating polygon domains.
* \ingroup	WlzPolyline
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>

/*!
* \return	Woolz error code.
* \ingroup	WlzPolyline
* \brief	Access to an polygon domain's integer vertices.
* \param	poly			Given polygon domain.
* \param	dstArySz		Destination pointer for the number of
* 					vertices.
* \param	dstPolyAry		Destination pointer for the vertices.
*/
WlzErrorNum	WlzPolyVertices2I(WlzPolygonDomain *poly,
				  int *dstArySz, WlzIVertex2 **dstPolyAry)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(poly == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstArySz == NULL) || (dstPolyAry == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstArySz = poly->nvertices;
    *dstPolyAry = poly->vtx;
  }
}

/*!
* \return	Woolz error code.
* \ingroup	WlzPolyline
* \brief	Access to an polygon domain's double precision vertices.
* \param	poly			Given polygon domain.
* \param	dstArySz		Destination pointer for the number of
* 					vertices.
* \param	dstPolyAry		Destination pointer for the vertices.
*/
WlzErrorNum	WlzPolyVertices2D(WlzPolygonDomain *poly,
				  int *dstArySz, WlzDVertex2 **dstPolyAry)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(poly == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((dstArySz == NULL) || (dstPolyAry == NULL))
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    *dstArySz = poly->nvertices;
    *dstPolyAry = (WlzDVertex2 *)(poly->vtx);
  }
}


