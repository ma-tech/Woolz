#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzTransform_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libWlz/WlzTransform.c
* \author       Bill Hill
* \date         April 2006
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
* \brief	Functions operating on Woolz transform unions.
* \ingroup	WlzTransform
*/

#include <Wlz.h>

/*!
* \return	Woolz error number.
* \ingroup	WlzTransform
* \brief	Free's the given Woolz transform.
* \param	tr			Given transform.
*/
WlzErrorNum	WlzFreeTransform(WlzTransform tr)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr.core == NULL)
  {
    errNum = WLZ_ERR_TRANSFORM_NULL;
  }
  else
  {
    switch(tr.core->type)
    {
      case WLZ_TRANSFORM_EMPTY:
        break;
      case WLZ_TRANSFORM_2D_AFFINE:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_REG:     /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_TRANS:   /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D_NOSHEAR: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_AFFINE:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_REG:     /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_TRANS:   /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_NOSHEAR:
	errNum = WlzFreeAffineTransform(tr.affine);
        break;
      case WLZ_TRANSFORM_2D_BASISFN:  /* FALLTHROUGH */
      case WLZ_TRANSFORM_2D5_BASISFN: /* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_BASISFN:
	errNum = WlzBasisFnFreeTransform(tr.basis);
        break;
      case WLZ_TRANSFORM_2D_MESH:	/* FALLTHROUGH */
	errNum = WlzMeshFreeTransform(tr.mesh);
        break;
      case WLZ_TRANSFORM_2D5_MESH:	/* FALLTHROUGH */
      case WLZ_TRANSFORM_3D_MESH:
	errNum = WLZ_ERR_UNIMPLEMENTED;
        break;
      default:
	errNum = WLZ_ERR_TRANSFORM_TYPE;
        break;
    }
  }
  return(errNum);
}
